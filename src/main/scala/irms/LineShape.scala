package irms {

    abstract class LineShape(expir:Array[Float], thirpeaks:Array[(Float,Float)]) {

        // TODO: What is the relationship between IR intensity in km/mol and TRANSMITTANCE?
        //       How do we convert these values?

        // information of single peak
        val nparams1:Int
        def f1(freq:Float,max:Float)(gparams:Array[Float],params1:Array[Float])(x:Float):Float
        def derivative1(freq:Float,max:Float)(gparams:Array[Float],params1:Array[Float])(x:Float):(Array[Float],Array[Float])

        // info of x
        // TODO: import values of xmin,xmax,xstep from 03_create_mid_struct_table.scala
        val xmin = 670 //included
        val xmax = 3702 //included
        val xstep = 4
        val xs = Range(xmin,xmax+xstep,xstep).toArray

        // Note: these functions are implemented assuming that the IR intensity and the
        //       TRANSMITTANCE have the same unit
        // TODO: fix this

        // theoretical spectrum vector
        // formula:
        // thirvec(xi) = sum( j from 0 to m, f1(freqj,maxj)(param1j)(xi) )
        def thirvec(gparams:Array[Float], params:Array[Float]):Array[Float] = {
            // split params into [param1]
            val param1s = params.sliding(nparams1,nparams1).toArray
            // apply (freq,max) to f1
            val fs = thirpeaks.map(j => f1(j._1,j._2) _)
            // apply (gparams,params1) to f1(freq,max)
            val fparams = fs.zip(param1s).map(j=>j._1(gparams,j._2))
            // calculate sum( j from 0 to m, f1(freq,max)(gparams,params1j)(xi) )
            xs.map( x => fparams.map(_(x)).reduce(_+_) )
        }

        // total loss:
        // formula:
        // euclidean_loss = sum( i from 0 to n, ( thir(xi)-expir(xi) ) ^ 2 )
        def euclidean_loss(gparams:Array[Float], params:Array[Float]):Float = {
            thirvec(gparams,params).zip(expir).map(j=>(j._1-j._2)*(j._1-j._2)).reduce(_+_)
        }

        // d(euclidean_loss)/d(parameters)
        // formula:
        // derivative_params  = sum( i from 0 to n, 2 * ( thir(xi)-expir(xi) ) * d(thir(xi))/d(parameters) )
        //                    = [ sum{ i from 0 to n, 2 * ( thir(xi)-expir(xi) ) * derivative1(freqj,maxj)(gparams,param1j)(xi)._2 } for j from 0 to m ]
        // derivative_gparams = sum( i from 0 to n, 2 * ( thir(xi)-expir(xi) ) * d(thir(xi))/d(gparameters) )
        //                    = sum{ i from 0 to n, j from 0 to m, 2 * ( thir(xi)-expir(xi) ) * derivative1(freqj,maxj)(gparams,param1j)(xi)._1 }
        def derivative_euclidean_loss(gparams:Array[Float], params:Array[Float]):(Array[Float],Array[Float]) = {
            // split params into [param1]
            val param1s = params.sliding(nparams1,nparams1).toArray
            // apply (freq,max) to derivative1
            val derivatives = thirpeaks.map(j => derivative1(j._1,j._2) _)
            // apply (gparams,params1) to derivative1(freq,max)
            val dparams = derivatives.zip(param1s).map(j=>j._1(gparams,j._2))
            // calculate 2 * ( thir(xi)-expir(xi) ) for all xi
            val diff2 = thirvec(gparams,params).zip(expir).map(j=>2*(j._1-j._2))
            // calculate elementwise plus/multiplies of two array:
            def cplus[C](a1:C,a2:C):C = ((a1,a2).asInstanceOf[(Any,Any)] match {
                case (a1:Array[Float],a2:Array[Float]) => a1.zip(a2).map(j=>j._1+j._2)
                case (a1:(Any,Any),a2:(Any,Any)) => (cplus[Any](a1._1,a2._1),cplus[Any](a1._2,a2._2))
            }).asInstanceOf[C]
            def cmult[C](collection:C,scalar:Float):C = (collection match {
                case c:Array[Float] => c.map(scalar*_)
                case c:(Any,Any) => (cmult[Any](c._1,scalar),cmult[Any](c._2,scalar))
            }).asInstanceOf[C]
            // calculate the sum{ i from 0 to n,.... }
            def sumi_calculator(d:Float=>(Array[Float],Array[Float])):(Array[Float],Array[Float]) = {
                val diff2deriv = diff2.zip(xs).map((j:(Float,Int))=>cmult(d(j._2),j._1))
                diff2deriv.reduce(cplus[(Array[Float],Array[Float])])
            }
            val sumi = dparams.map(sumi_calculator)
            ( sumi.map(_._2).reduce(cplus[Array[Float]]), sumi.flatMap(_._1) )
        }

    }

    // private class ForcedOscillatior(expir:Array[Float], thir:Array[(Float,Float)]) extends LineShape {
    //     override val nparams1:Int = 1
    //
    //     // formula A = F / sqrt( (f0^2-f^2)^2 + 4*beta^2*f^2 )
    //     //         E = A^2
    //     def f1(freq:Float,max:Foat)(gparams:Array[Float], params1:Array[Float])(x:Float):Float = {
    //         val beta = params1(0)
    //         val b2f24 = 4 * beta*beta * freq*freq
    //         val FF = max * b2f24
    //         val diff2freq = freq*freq - x*x
    //         FF / ( diff2freq*diff2freq + b2f24 )
    //     }
    //
    //     // formula dE/d(beta) = ...
    //     def derivative1(freq:Float,max:Float)(gparams:Array[Float], params1:Array[Float])(x:Float):Array[Float] = {
    //         //TODO: implement
    //     }
    // }
}
