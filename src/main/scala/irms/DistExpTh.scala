import org.apache.spark.sql._
import org.apache.spark._
import org.apache.spark.sql.functions._
import scala.math._

package irms {

    object DistExpTh {

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

        case class MyExpIR(smiles:String,vec:Array[Double])
        case class MyThIR(smiles:String, freqs:Array[(Double,Double)])
        case class MyDis(th_smiles:String,exp_smiles:String,expvec:Array[Double],thvec:Array[Double],distance:Double)

        def main(args: Array[String]):Unit = {
            // read tables
            val session = SparkSession.builder.appName("d(thir,expir)").getOrCreate()
            import session.implicits._
            val path = "/home/gaoxiang/irms/create-dataset-for-ir/outputs/tables"
            val mid_structure = session.read.parquet(path+"/mid_structure").as[MIDStruct]
            val expir = session.read.parquet(path+"/expir").as[ExpIRAndState]
            val thir = session.read.parquet(path+"/thir").as[TheoreticalIR]
            // get all gas phases
            val thir_b3lyp631gd = thir.filter(_.method=="B3LYP/6-31G*")
            val expir_selected = expir.filter(_.state=="gas")
                                      .joinWith(mid_structure,mid_structure("mid")===expir("mid"))
                                      .map(j=>MyExpIR(j._2.smiles,j._1.vec)).limit(4)
            val thir_selected = thir_b3lyp631gd.joinWith(expir_selected,expir_selected("smiles")===thir_b3lyp631gd("smiles"))
                                               .map(j=>MyThIR(j._1.smiles,j._1.freqs)).limit(4)
            // calculate distances
            def calculate_distance(j:(MyThIR,MyExpIR)):MyDis = {
                val (th,exp) = j
                val calculator = new Lorentz(th.freqs.toSeq)
                val d = calculator.loss(exp.vec.toSeq)
                val thvec = LineShapeHelpers.vec(calculator.thir(calculator.final_gparams.get,calculator.final_params.get))
                MyDis(th.smiles,exp.smiles,exp.vec,thvec.toArray,d)
            }
            val distances = thir_selected.joinWith(expir_selected,expr("1>0")).repartition(400).map(calculate_distance _)
            distances.show()

            // write result to file
            val resultpath = "/home/gaoxiang/irms/exp-th-distance"
            distances.write.parquet(resultpath+"/distances")
        }
    }
}
