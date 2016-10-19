package irms {

	object DistExpTh {
		import scala.collection.immutable.{SortedSet,SortedMap}
		import scala.math.{abs,sqrt}

		// set up breeze
		import scala.language.implicitConversions
		import breeze.optimize._
		import breeze.optimize.StochasticGradientDescent._
		import breeze.linalg._
		type DV = DenseVector[Double]
		type DF = DiffFunction[DV]
		implicit def seq2densevec(s:Seq[Double]):DV = new DV(s.toArray)
		implicit def densevec2seq(v:DV):Seq[Double] = v.data.toSeq

		val width = 7.0
		val uniform_scaling = 0.96
		private class MyUniform(peaks:Seq[(Double,Double)]) extends LineShape(peaks) with UniformScaled with CosineLoss with ZeroBaseline with Lorentzian {
			val scaling = uniform_scaling
			val w = width
		}
		private class MyPiecewise(s:SortedMap[Double,Double], val smiles:String,peaks:Seq[(Double,Double)]) extends LineShape(peaks) with PiecewiseScaled with CosineLoss with ZeroBaseline with Lorentzian {
			val scalings = s
			val w = width
		}

		def dlorentz(freq:Double,intensity:Double)(x:Double):Double = {
			val d = x - freq
			val den1 = d * d + width * width
			- 2 * d * intensity * width / (den1*den1)
		}

		val uppers:SortedSet[Double] = SortedSet(1000,1500,2000,2500,3000,3500,4000)
		val init_scalings:Seq[Double] = uppers.toSeq.map(j=>uniform_scaling)

		def train_scalings1(peaks:Array[(Double,Double)],expir:Array[Double]):Array[Double] = {
			val difffunc = new DF {
				val upper_peaks_map = SortedMap(SortedSet(peaks: _*).groupBy(j=>uppers.from(j._1).head).toArray: _*)
				val empty_uppers = (uppers &~ upper_peaks_map.keys.toSet).map((_,SortedSet[(Double,Double)]())).toMap
				val peaks_by_range:Seq[SortedSet[(Double,Double)]] = (upper_peaks_map++empty_uppers).values.toSeq

				def calculate(scalings:DV):(Double,DV) = {
					val scaling_map = SortedMap(uppers.toSeq.zip(scalings).toArray: _*)
					val lineshape = new MyPiecewise(scaling_map,"",peaks)
					val thir = lineshape.thir
					val alpha = sqrt(thir.map(a=>a*a).sum)
					val beta = (thir,expir).zipped.map(_*_).sum
					val gamma = sqrt(expir.map(a=>a*a).sum)
					val uniform_lineshape = new MyUniform(peaks)
					val dthir_ds = (scalings,peaks_by_range).zipped.map( (s,peaks)=>
						Params.X.xs.map( x =>
							-peaks.unzip.zipped.map((f:Double,a:Double)=>f*dlorentz(s*f,a)(x)).sum
						)
					)
					val dalpha_ds = dthir_ds.map(dthir=>(dthir,thir).zipped.map(_*_).sum/alpha)
					val dbeta_ds = dthir_ds.map(dthir=>(dthir,expir).zipped.map(_*_).sum)
					val ds = (dalpha_ds,dbeta_ds).zipped.map((dalpha,dbeta) => -(alpha*dbeta-beta*dalpha)/(gamma*alpha*alpha))
					(lineshape.loss(expir)-uniform_lineshape.loss(expir),ds)
				}
			}

			// val optimizer = new SimpleSGD[DV](0.001)
			val optimizer = new LBFGS[DV]()
			val optimum = optimizer.minimize(difffunc,init_scalings)
			optimum.toArray
		}

		case class MyThExp(mid:String, smiles:String, vec:Array[Double], freqs:Array[(Double,Double)])
		case class MyExp(smiles:String, vec:Array[Double])
		case class MyTh(smiles:String, freqs:Array[(Double,Double)])
		case class MyThExpVec(mid:String, smiles:String, expvec:Array[Double], thvec:Array[Double])
		case class MyRD(exp_smiles:String,rank:Int,distance:Double,mistakes:Array[String],distances:Array[Double])

		def main(args: Array[String]):Unit = {

			import org.apache.spark.sql._
			import org.apache.spark._
			import org.apache.spark.sql.functions._
			import Env.spark.implicits._

			// pick data
			val mid_structure = TableManager.getOrCreate(MIDStruct)
			val expir = TableManager.getOrCreate(ExpIRAndState).filter(_.state=="gas")
			val thir = TableManager.getOrCreate(TheoreticalIR).filter(_.method=="B3LYP/6-31G*")
			val thexp = expir.join(mid_structure,"mid").join(thir,"smiles")

			// train scaling factors
			val splited_thexp = thexp.randomSplit(Array(0.02,0.98))
			val optimized_scalings = splited_thexp(0).as[MyThExp].map(j=>train_scalings1(j.freqs,j.vec))
			println("optimized scalings:")
			optimized_scalings.show()
			val scaling_count = optimized_scalings.count()
			def sum_scaling(a:Array[Double],b:Array[Double]):Array[Double] = {
				(a.toSeq,b.toSeq).zipped.map(_+_).toArray
			}
			val avg_scalings = optimized_scalings.reduce(sum_scaling _).toSeq.map(_/scaling_count)
			val scaling_map = SortedMap(uppers.zip(avg_scalings).toArray:_*)

			// output vectorized thir
			val thexpvec = splited_thexp(1).as[MyThExp].map(j=>MyThExpVec(j.mid,j.smiles,j.vec,new MyPiecewise(scaling_map,j.smiles,j.freqs).thir.toArray))
			thexpvec.show()
			thexpvec.write.parquet(Env.tables+"/thexpvec")

			// calculate distances
			val thir_local = splited_thexp(1).select("smiles","freqs").dropDuplicates("smiles").as[MyTh].rdd.map(j=>new MyPiecewise(scaling_map,j.smiles,j.freqs)).collect().toSeq
			def expir2rd(exp:MyExp):MyRD = {
				val distance_smiles = thir_local.map(j=>(j.loss(exp.vec.toSeq),j.smiles)).sorted
				val ((distance,_),rank) = distance_smiles.zipWithIndex.find(j=>j._1._2==exp.smiles).get
				val mistakes = distance_smiles.splitAt(rank)._1.unzip._2.toArray
				MyRD(exp.smiles,rank,distance,mistakes,distance_smiles.unzip._1.toArray)
			}
			val myexp = splited_thexp(1).select("smiles","vec").dropDuplicates("smiles").as[MyExp]
			val count = thir_local.length
			val rd = myexp.repartition(count/20).map(expir2rd _)
			rd.show()
			println("total number of data: "+count)
			rd.write.parquet(Env.tables+"/distances")
		}
	}
}
