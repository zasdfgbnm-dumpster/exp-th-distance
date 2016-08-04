import org.apache.spark.sql._
import org.apache.spark._
import scala.math._

package irms {

    object Tests {
        private class TestGaussian(thirpeaks:Array[(Float,Float)]) extends LineShape(Array[Float](),thirpeaks) {
            override val nparams1:Int = 0
            val s = 20
            def baseline(gparams:Array[Float])(x:Float):Float = x * 0.01f
            def dbaseline(gparams:Array[Float])(x:Float):Array[Float] = gparams
            def f1(freq:Float,max:Float)(gparams:Array[Float], params1:Array[Float])(x:Float):Float = {
                val xx = (x-freq)/s
                max * exp(-xx*xx).toFloat
            }
            def df1(freq:Float,max:Float)(gparams:Array[Float], params1:Array[Float])(x:Float):(Array[Float],Array[Float]) = (gparams,gparams)
        }

        def main(args: Array[String]):Unit = {
            val session = SparkSession.builder.appName("draw_thir").getOrCreate()
            import session.implicits._
            val path = "/home/gaoxiang/irms/create-dataset-for-ir/outputs/06/thir"

            // draw thir
            val thir = session.read.parquet(path).as[TheoreticalIR]
            val thirvec = thir.map(j=>new TestGaussian(j.freqs).thirvec(Array[Float](),Array[Float]()))
            val thirstr = thirvec.map(a=>a.map(_.toString).reduce(_+" "+_))
            val thirstrs = thirstr.reduce(_+"\n"+_)
            println(thirstrs)
        }
    }
}
