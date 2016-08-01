import org.apache.spark.sql._
import org.apache.spark._

package irms {

    object DistExpTh {

        private case class ParameterError(msg:String) extends Throwable

        private abstract class LineShape(expir:Array[Float], thir:Array[Float]) {

            // TODO: What is the relationship between IR intensity in km/mol and TRANSMITTANCE?
            //       How do we convert these values?

            // information of single peak
            val nparams1:Int
            def f1(freq:Float,max:Foat)(params1:Array[Float])(x:Float):Float
            def derivative1(freq:Float,max:Foat)(params1:Array[Float])(x:Float):Array[Float]

            // info of x
            // TODO: import values of xmin,xmax,xstep from 03_create_mid_struct_table.scala
            val xs = Range(xmin,xmax+xstep,xstep)

            // Note: these functions are implemented assuming that the IR intensity and the
            //       TRANSMITTANCE have the same unit
            // TODO: fix this

            // theoretical spectrum vector
            // formula:
            // thir(xi) = sum( j from 0 to m, f1(freqj,maxj)(xi)(param1j) )
            def thir(params:Array[Float]):Array[Float] = {
                // split params into [param1]
                val param1s = params.sliding(nparams1,nparams1)
                val fs = thir.map(j => f1(j._1,j._2))
                val fparams = fs.zip(param1s).map(j=>j._1(j._2))
                xs.map( x => fparams.map(_(x)).reduce(_+_) )
            }

            // total loss:
            // formula:
            // euclidean_loss = sum( i from 0 to n, ( thir(xi)-expir(xi) ) ^ 2 )
            def euclidean_loss(params:Array[Float]):Float = {
                thir(params).zip(expir).map(j=>(j._1-j._2)*(j._1-j._2)).reduce(_+_)
            }

            // d(euclidean_loss)/d(parameters)
            // formula:
            // derivative = sum( i from 0 to n, 2 * ( thir(xi)-expir(xi) ) * d(thir(xi))/d(parameters) )
            //            = [ sum{ i from 0 to n, 2 * ( thir(xi)-expir(xi) ) * derivative1(freqj,maxj)(param1j)(xi) } for j from 0 to m ]
            //
            def derivative_euclidean_loss(params:Array[Float]):Array[Float] = {
                val param1s = params.sliding(nparams1,nparams1)
                val derivatives = thir.map(j => derivative1(j._1,j._2))
                val dparams = fs.zip(param1s).map(j=>j._1(j._2))
                val diff2 = thir(params).zip(expir).map(j=>2*(j._1-j._2))
                def sum(d:Float=>Array[Float]):Array[Float] = {
                    val diff2deriv = diff2.zip(xs).map(j=>d(j._2).map(j._1*_))
                    diff2deriv.reduce(_.zip(_).map(a=>a._1+a._2))
                }
                dparams.flatMap(sum)
            }

        }

        private class ForcedOscillatior(expir:Array[Float], thir:Array[Float]) extends LineShape {
            override val nparams1:Int = 1

            // formula A = F / sqrt( (f0^2-f^2)^2 + 4*beta^2*f^2 )
            //         E = A^2
            override def f1(freq:Float,max:Foat)(params1:Array[Float])(x:Float):Float = {
                val beta = params(0)
                val b2f24 = 4 * beta*beta * freq*freq
                val FF = max * b2f24
                val diff2freq = freq*freq - x*x
                FF / ( diff2freq*diff2freq + b2f24 )
            }

            // formula dE/d(beta) = ...
            override def derivative1(freq:Float,max:Foat)(params1:Array[Float])(x:Float):Array[Float] = {
                //TODO: implement
            }
        }

        def main(args: Array[String]):Unit = {
            val session = SparkSession.builder.appName("03_create_mid_struct_table").getOrCreate()
            import session.implicits._
            val path = "/home/gaoxiang/create-dataset-for-ir/outputs/"
        }
    }
}
