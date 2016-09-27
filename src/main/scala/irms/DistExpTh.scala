package irms {

    object DistExpTh {

        private class Gaussian(peaks:Seq[(Double,Double)]) extends LineShape(peaks) with Optimizer.SimpleSGDOptimizer with Loss.Euclidean.Optmized with Functions.ZeroBaseline.Optimized {
            import scala.math._
            val r = scala.util.Random
            val nparams1:Int = 1
            def gparams:Seq[Double] = Seq(1.0/peaks.unzip._2.max*(0.7+1.3*r.nextDouble))
            def params:Seq[Double] = Range(0,peaks.length*nparams1).map(j=>r.nextDouble*10+5)
            def f1(freq:Double,intensity:Double,gparams:Seq[Double],params1:Seq[Double])(x:Double):Double = {
                val scaled_freq = freq * 0.960
                val a = intensity * gparams(0)
                val w = params1(0)
                val d = ( x - scaled_freq )
                val scaled_d = d/w
                (a/w) * exp(-scaled_d*scaled_d)
            }
            def df1(freq:Double,intensity:Double,gparams:Seq[Double],params1:Seq[Double])(x:Double):(Seq[Double],Seq[Double]) = {
                val f1val = f1(freq,intensity,gparams,params1)(x)
                val scaled_freq = freq * 0.960
                val a = intensity * gparams(0)
                val w = params1(0)
                val d = ( x - scaled_freq )
                val df1_da = f1val/a
                val df1_dw = f1val*(2*d*d-w*w)/(w*w*w)
                (Seq(df1_da),Seq(df1_dw))
            }
        }

        private class Lorentz(peaks:Seq[(Double,Double)]) extends LineShape(peaks) with Optimizer.SimpleSGDOptimizer with Loss.Euclidean.Optmized with Functions.ZeroBaseline.Optimized {
            val r = scala.util.Random
            val nparams1:Int = 1
            def gparams:Seq[Double] = Seq(1.0/peaks.unzip._2.max*(0.7+1.3*r.nextDouble))
            def params:Seq[Double] = Range(0,peaks.length*nparams1).map(j=>r.nextDouble*20+5)
            def f1(freq:Double,intensity:Double,gparams:Seq[Double],params1:Seq[Double])(x:Double):Double = {
                val scaled_freq = freq * 0.960
                val a = intensity * gparams(0)
                val w = params1(0)
                val d = ( x - scaled_freq )
                a * w / ( d * d + w * w )
            }
            def df1(freq:Double,intensity:Double,gparams:Seq[Double],params1:Seq[Double])(x:Double):(Seq[Double],Seq[Double]) = {
                val f1val = f1(freq,intensity,gparams,params1)(x)
                val scaled_freq = freq * 0.960
                val a = intensity * gparams(0)
                val w = params1(0)
                val d = ( x - scaled_freq )
                val df1_da = f1val/a
                val df1_dw = (f1val*f1val/a)*(d*d/(w*w)-1)
                (Seq(df1_da),Seq(df1_dw))
            }
        }

        val path = "/ufrc/roitberg/qasdfgtyuiop/tables"
        val resultpath = "/ufrc/roitberg/qasdfgtyuiop/exp-th-dis"
        // val path = "tables"
        // val resultpath = "."

        case class MyExpIR(smiles:String,vec:Array[Double])
        case class MyThIR(smiles:String, freqs:Array[(Double,Double)])
        case class MyRD(exp_smiles:String,rank:Int,distance:Double,mistakes:Array[String],distances:Array[Double])

        def main(args: Array[String]):Unit = {

            import org.apache.spark.sql._
            import org.apache.spark._
            import org.apache.spark.sql.functions._

            val session = SparkSession.builder.appName("d(thir,expir)").getOrCreate()
            import session.implicits._

            // pick data
            val mid_structure = session.read.parquet(path+"/mid_structure").as[MIDStruct]
            val expir = session.read.parquet(path+"/expir").as[ExpIRAndState]
            val thir = session.read.parquet(path+"/thir").as[TheoreticalIR]

            // select data
            def normalize_vec(expir:MyExpIR):MyExpIR = {
                MyExpIR(expir.smiles,LineShapeHelpers.cmult(expir.vec.toSeq,1.0/expir.vec.max).toArray)
            }
            val thir_b3lyp631gd = thir.filter(_.method=="B3LYP/6-31G*")
            val expir_gas = expir.filter(_.state=="gas")
            val thir_selected = thir_b3lyp631gd.map(j=>MyThIR(j.smiles,j.freqs))
            val thir_local = thir_selected.collect().toSeq
            val expir_converted = expir_gas.joinWith(mid_structure,mid_structure("mid")===expir("mid"))
                                           .map(j=>MyExpIR(j._2.smiles,j._1.vec))
            val expir_normalized = expir_converted.map(normalize_vec _)
            val expir_included = expir_normalized.joinWith(thir_selected,expir_normalized("smiles")===thir_selected("smiles")).map(_._1)
            val count_expir = expir_included.count().toInt

            // calculate distances
            def dis(exp:MyExpIR,th:MyThIR):Double = {
                new Gaussian(th.freqs.toSeq).loss(exp.vec.toSeq)
            }
            val r = scala.util.Random
            def dis_debug(exp:MyExpIR,th:MyThIR):Double = {
                scala.math.abs(exp.smiles.length-th.smiles.length) + r.nextDouble * 3
            }
            def expir2rd(exp:MyExpIR):MyRD = {
                val distance_smiles = thir_local.map(j=>(dis_debug(exp,j),j.smiles)).sorted
                val ((distance,_),rank) = distance_smiles.zipWithIndex.find(j=>j._1._2==exp.smiles).get
                val mistakes = distance_smiles.splitAt(rank)._1.unzip._2.toArray
                MyRD(exp.smiles,rank,distance,mistakes,distance_smiles.unzip._1.toArray)
            }
            val rd = expir_included.repartition(count_expir).map(expir2rd _)
            rd.show()
            rd.write.parquet(resultpath+"/rd")
        }
    }
}
