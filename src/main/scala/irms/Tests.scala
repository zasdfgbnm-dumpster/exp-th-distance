import org.apache.spark.sql._
import org.apache.spark._
import scala.math._

package irms {

    object Tests {

        import X.xs
        import LineShapeHelpers.vec

        private class TestGaussian(peaks:Seq[(Float,Float)]) extends LineShape(peaks) {
            override val nparams1:Int = 0
            val s = 20
            def baseline(gparams:Seq[Float])(x:Float):Float = x * 0.01f
            def f1(freq:Float,max:Float,gparams:Seq[Float], params1:Seq[Float])(x:Float):Float = {
                val xx = (x-freq)/s
                max * exp(-xx*xx).toFloat
            }
        }

        def main(args: Array[String]):Unit = {
            val session = SparkSession.builder.appName("draw_thir").getOrCreate()
            import session.implicits._
            val path = "/home/gaoxiang/irms/create-dataset-for-ir/outputs/06/thir"

            // draw thir
            val thir = session.read.parquet(path).as[TheoreticalIR]
            val thirvec = thir.map(j=>vec(new TestGaussian(j.freqs).thir(Seq[Float](),Seq[Float]())).toArray)
            val thirstr = thirvec.map(a=>a.map(_.toString).reduce(_+" "+_))
            val thirstrs = thirstr.take(20).reduce(_+"\n"+_)
            println(thirstrs)
        }
    }
}
