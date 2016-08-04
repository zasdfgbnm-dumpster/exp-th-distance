package irms {

    abstract class LineShape(expir:Array[Float], thirpeaks:Array[(Float,Float)]) {

        // information of single peak
        // nparams1: number of local parameters that f1 and derivative1 take
        val nparams1:Int
        // f1: the line shape of a single peak, gparams are global parameters, params1 are local parameters
        def f1(freq:Float,max:Float)(gparams:Array[Float],params1:Array[Float])(x:Float):Float
        // derivative1: derivative of f1 with gparams and params1.
        // return_value._1 are derivative with gparams, return_value._2 are with params1
        def derivative1(freq:Float,max:Float)(gparams:Array[Float],params1:Array[Float])(x:Float):(Array[Float],Array[Float])

        val xs = Range(X.xmin,X.xmax+X.xstep,X.xstep).toArray

        // theoretical spectrum vector
        // formula:
        // thirvec(xi) = sum( j from 0 to m, f1(freqj,maxj)(param1j)(xi) )
        def thirvec(gparams:Array[Float], params:Array[Float]):Array[Float] = {
            // split params into [param1]
            val param1s = if(nparams1==0) {
                thirpeaks.map(j=>Array[Float]())
            } else params.sliding(nparams1,nparams1).toArray
            // apply (freq,max) to f1
            val fs = thirpeaks.map(j => f1(j._1,j._2) _)
            // apply (gparams,params1) to f1(freq,max)
            val fparams = fs.zip(param1s).map(j=>j._1(gparams,j._2))
            // calculate sum( j from 0 to m, f1(freq,max)(gparams,params1j)(xi) )
            xs.map( x => fparams.map(_(x)).reduce(_+_) )
        }

        object loss {
            // total loss:
            // formula:
            // euclidean_loss = sum( i from 0 to n, ( thir(xi)-expir(xi) ) ^ 2 )
            def euclidean(gparams:Array[Float], params:Array[Float]):Float = {
                thirvec(gparams,params).zip(expir).map(j=>(j._1-j._2)*(j._1-j._2)).reduce(_+_)
            }

            // d(euclidean)/d(parameters)
            // formula:
            // derivative_gparams = sum( i from 0 to n, 2 * ( thir(xi)-expir(xi) ) * d(thir(xi))/d(gparameters) )
            //                    = sum{ i from 0 to n, j from 0 to m, 2 * ( thir(xi)-expir(xi) ) * derivative1(freqj,maxj)(gparams,param1j)(xi)._1 }
            // derivative_params  = sum( i from 0 to n, 2 * ( thir(xi)-expir(xi) ) * d(thir(xi))/d(parameters) )
            //                    = [ sum{ i from 0 to n, 2 * ( thir(xi)-expir(xi) ) * derivative1(freqj,maxj)(gparams,param1j)(xi)._2 } for j from 0 to m ]
            def deuclidean(gparams:Array[Float], params:Array[Float]):(Array[Float],Array[Float]) = {
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
                ( sumi.map(_._1).reduce(cplus[Array[Float]]), sumi.flatMap(_._2) )
            }
        }

    }

}
