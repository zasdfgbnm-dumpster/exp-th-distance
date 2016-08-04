package irms {

    abstract class LineShape(expir:Array[Float], thirpeaks:Array[(Float,Float)]) {

        // information of single peak

        // nparams1: number of local parameters that f1 and derivativef1 take
        val nparams1:Int

        // baseline: baseline of the spectrum
        // gparams are global parameters shared with f1
        def baseline(gparams:Array[Float])(x:Float):Float
        def derivativeb(gparams:Array[Float])(x:Float):Array[Float]

        // f1: the line shape of a single peak,
        // gparams are global parameters (shared with baseline), params1 are local parameters
        def f1(freq:Float,max:Float)(gparams:Array[Float],params1:Array[Float])(x:Float):Float
        // derivativef1: derivative of f1 with gparams and params1.
        // return_value._1 are derivative with gparams, return_value._2 are with params1
        def derivativef1(freq:Float,max:Float)(gparams:Array[Float],params1:Array[Float])(x:Float):(Array[Float],Array[Float])

        val xs = Range(X.xmin,X.xmax+X.xstep,X.xstep).toArray.map(1.0f*_)

        // calculate elementwise plus/multiplies of two array/tuple2:
        private def cplus[C](a1:C,a2:C):C = ((a1,a2).asInstanceOf[(Any,Any)] match {
            case (a1:Array[Float],a2:Array[Float]) => a1.zip(a2).map(j=>j._1+j._2)
            case (a1:(Any,Any),a2:(Any,Any)) => (cplus[Any](a1._1,a2._1),cplus[Any](a1._2,a2._2))
        }).asInstanceOf[C]
        private def cmult[C](collection:C,scalar:Float):C = (collection match {
            case c:Array[Float] => c.map(scalar*_)
            case c:(Any,Any) => (cmult[Any](c._1,scalar),cmult[Any](c._2,scalar))
        }).asInstanceOf[C]

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
            val peaks = xs.map( x => fparams.map(_(x)).reduce(_+_) )
            // calculate baseline(gparams)(xi)
            val blvec = xs.map( baseline(gparams)_ )
            // calculate baseline + peaks
            cplus(peaks,blvec)
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
            //                    = sum{ i from 0 to n, 2 * ( thir(xi)-expir(xi) ) * ( derivativeb + sum { j from 0 to m, derivativef1(freqj,maxj)(gparams,param1j)(xi)._1 })}
            // derivative_params  = sum( i from 0 to n, 2 * ( thir(xi)-expir(xi) ) * d(thir(xi))/d(parameters) )
            //                    = [ sum{ i from 0 to n, 2 * ( thir(xi)-expir(xi) ) * derivativef1(freqj,maxj)(gparams,param1j)(xi)._2 } for j from 0 to m ]
            def deuclidean(gparams:Array[Float], params:Array[Float]):(Array[Float],Array[Float]) = {
                // split params into [param1]
                val param1s = params.sliding(nparams1,nparams1).toArray
                // apply (freq,max) to derivativef1
                val derivatives = thirpeaks.map(j => derivativef1(j._1,j._2) _)
                // apply (gparams,params1) to derivativef1(freq,max)
                val dparams = derivatives.zip(param1s).map(j=>j._1(gparams,j._2))
                // calculate 2 * ( thir(xi)-expir(xi) ) for all xi
                val diff2 = thirvec(gparams,params).zip(expir).map(j=>2*(j._1-j._2))
                // calculate the sum{ i from 0 to n,.... } for peaks
                def sumi_calculator(d:Float=>(Array[Float],Array[Float])):(Array[Float],Array[Float]) = {
                    val diff2deriv = diff2.zip(xs).map(j=>cmult(d(j._2),j._1))
                    diff2deriv.reduce(cplus[(Array[Float],Array[Float])])
                }
                val sumi = dparams.map(sumi_calculator)
                val dgparams_peaks = sumi.map(_._1).reduce(cplus[Array[Float]])
                val dparams_peaks = sumi.flatMap(_._2)
                // calculate the sum{ i from 0 to n,.... } for baseline
                val dgparams_baseline = diff2.zip(xs).map(j=>cmult(baseline(gparams)(j._2),j._1))
                // final result
                ( cplus(dgparams_peaks,dgparams_baseline) , dparams_peaks )
            }
        }

    }

}
