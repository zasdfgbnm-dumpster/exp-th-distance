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

        case class MyExpIR(smiles:String,vec:Array[Double])
        case class MyThIR(smiles:String, freqs:Array[(Double,Double)])
        case class MyDis(th_smiles:String,exp_smiles:String,expvec:Array[Double],thvec:Array[Double],distance:Double)
        case class MyRD(th_smiles:String,rank:Int,d:Double,mistakes:Array[String],ds:Array[Double])

        def main(args: Array[String]):Unit = {
            import org.apache.spark.sql._
            import org.apache.spark._
            import org.apache.spark.sql.functions._

            // read tables
            val session = SparkSession.builder.appName("d(thir,expir)").getOrCreate()
            import session.implicits._
            val mid_structure = session.read.parquet(path+"/mid_structure").as[MIDStruct]
            val expir = session.read.parquet(path+"/expir").as[ExpIRAndState]
            val thir = session.read.parquet(path+"/thir").as[TheoreticalIR]

            // select data
            def normalize_vec(expir:MyExpIR):MyExpIR = {
                MyExpIR(expir.smiles,LineShapeHelpers.cmult(expir.vec.toSeq,1.0/expir.vec.max).toArray)
            }
            val thir_b3lyp631gd = thir.filter(_.method=="B3LYP/6-31G*")
            val expir_selected = expir.filter(_.state=="gas")
                                      .joinWith(mid_structure,mid_structure("mid")===expir("mid"))
                                      .map(j=>MyExpIR(j._2.smiles,j._1.vec))
                                      .map(normalize_vec _)
            val thir_selected = thir_b3lyp631gd.joinWith(expir_selected,expir_selected("smiles")===thir_b3lyp631gd("smiles"))
                                               .map(j=>MyThIR(j._1.smiles,j._1.freqs))

            // calculate distances
            def calculate_distance(j:(MyThIR,MyExpIR)):MyDis = {
                val (th,exp) = j
                val calculator = new Gaussian(th.freqs.toSeq)
                val d = calculator.loss(exp.vec.toSeq)
                val thvec = LineShapeHelpers.vec(calculator.thir(calculator.final_gparams.get,calculator.final_params.get))
                MyDis(th.smiles,exp.smiles,exp.vec,thvec.toArray,d)
            }
            val distances = thir_selected.joinWith(expir_selected,expr("true")).repartition(4000).map(calculate_distance _)
            distances.show()
            distances.write.parquet(resultpath+"/distances")

            // plot distance vs rank
            val ed = distances.rdd.map(j=>(j.th_smiles,(j.distance,j.exp_smiles))).groupByKey()
            def ed2rd(ed:(String,Iterable[(Double,String)])):MyRD = {
                val ed2s = ed._2.toSeq
                val th_smiles = ed._1
                val l = ed2s.sorted
                val ((d,_),r) = l.zipWithIndex.find(j=>j._1._2==th_smiles).get
                val mistakes = ed2s.splitAt(r)._1.unzip._2.toArray
                MyRD(th_smiles,r,d,mistakes,l.unzip._1.toArray)
            }
            val rd = ed.map(ed2rd).toDS
            rd.show()
            rd.write.parquet(resultpath+"/rd")
        }
    }
}
