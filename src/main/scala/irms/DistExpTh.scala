package irms {

	object DistExpTh {

		private class MyLineshape(val smiles:String,peaks:Seq[(Double,Double)]) extends LineShape(peaks) with UniformScaled with CosineLoss with ZeroBaseline with Lorentzian {
			val scaling = 0.96
			val w = 7.0
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
			val spark = Env.spark
			import spark.implicits._

			// pick data
			val mid_structure = TableManager.getOrCreate(MIDStruct).as[MIDStruct]
			val expir = TableManager.getOrCreate(ExpIRAndState).as[ExpIRAndState].filter(_.state=="gas")
			val thir = TableManager.getOrCreate(TheoreticalIR).as[TheoreticalIR].filter(_.method=="B3LYP/6-31G*")
			val thexp = expir.join(mid_structure,"mid").join(thir,"smiles")

			// output vectorized thir
			val thexpvec = thexp.as[MyThExp].map(j=>MyThExpVec(j.mid,j.smiles,j.vec,new MyLineshape(j.smiles,j.freqs).thir.toArray))
			thexpvec.show()
			thexpvec.write.parquet(Env.tables+"/thexpvec")

			// calculate distances
			val thir_local = thexp.select("smiles","freqs").dropDuplicates("smiles").as[MyTh].rdd.map(j=>new MyLineshape(j.smiles,j.freqs)).collect().toSeq
			def expir2rd(exp:MyExp):MyRD = {
				val distance_smiles = thir_local.map(j=>(j.loss(exp.vec.toSeq),j.smiles)).sorted
				val ((distance,_),rank) = distance_smiles.zipWithIndex.find(j=>j._1._2==exp.smiles).get
				val mistakes = distance_smiles.splitAt(rank)._1.unzip._2.toArray
				MyRD(exp.smiles,rank,distance,mistakes,distance_smiles.unzip._1.toArray)
			}
			val myexp = thexp.select("smiles","vec").dropDuplicates("smiles").as[MyExp]
			val count = thir_local.length
			val rd = myexp.repartition(count).map(expir2rd _)
			rd.show()
			println("total number of data: "+count)
			rd.write.parquet(Env.tables+"/distances")
		}
	}
}
