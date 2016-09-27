package irms {

	object DistExpTh {

		private class MyLineshape(peaks:Seq[(Double,Double)]) extends LineShape(peaks) with UniformScaled with CosineLoss with ZeroBaseline with Lorentzian {
			val scaling = 0.96
			val w = 7.0
		}

		case class MyExpIR(smiles:String,vec:Array[Double])
		case class MyThIR(smiles:String, freqs:Array[(Double,Double)])
		case class MyRD(exp_smiles:String,rank:Int,distance:Double,mistakes:Array[String],distances:Array[Double])

		def main(args: Array[String]):Unit = {

			import org.apache.spark.sql._
			import org.apache.spark._
			import org.apache.spark.sql.functions._
			import Env.spark.implicits._

			// pick data
			val mid_structure = TableManager.getOrCreate(MIDStruct)
			val expir = TableManager.getOrCreate(ExpIRAndState)
			val thir = TableManager.getOrCreate(TheoreticalIR)
			val thir_b3lyp631gd = thir.filter(_.method=="B3LYP/6-31G*")
			val expir_gas = expir.filter(_.state=="gas")
			val thir_selected = thir_b3lyp631gd.map(j=>MyThIR(j.smiles,j.freqs))
			val thir_local = thir_selected.collect().toSeq
			val expir_converted = expir_gas.joinWith(mid_structure,mid_structure("mid")===expir("mid"))
										   .map(j=>MyExpIR(j._2.smiles,j._1.vec))
			val expir_included = expir_converted.joinWith(thir_selected,expir_converted("smiles")===thir_selected("smiles")).map(_._1)
			val count_expir = expir_included.count().toInt

			// calculate distances
			def dis(exp:MyExpIR,th:MyThIR):Double = {
				new MyLineshape(th.freqs.toSeq).loss(exp.vec.toSeq)
			}
			val r = scala.util.Random
			def expir2rd(exp:MyExpIR):MyRD = {
				val distance_smiles = thir_local.map(j=>(dis(exp,j),j.smiles)).sorted
				val ((distance,_),rank) = distance_smiles.zipWithIndex.find(j=>j._1._2==exp.smiles).get
				val mistakes = distance_smiles.splitAt(rank)._1.unzip._2.toArray
				MyRD(exp.smiles,rank,distance,mistakes,distance_smiles.unzip._1.toArray)
			}
			val rd = expir_included.repartition(count_expir).map(expir2rd _)
			rd.show()
			//rd.write.parquet(Env.tables+"/rd")
		}
	}
}
