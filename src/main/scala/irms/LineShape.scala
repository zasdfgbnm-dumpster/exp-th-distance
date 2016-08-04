package irms {

    object LineShapeHelpers {

        // x values
        val xs = Range(X.xmin,X.xmax+X.xstep,X.xstep).toArray.map(1.0f*_)

        // calculate elementwise plus/multiplies of two array/tuple2:
        // only tuple2 and array of Float are supported now...
        def cplus[C](a1:C,a2:C):C = ((a1,a2).asInstanceOf[(Any,Any)] match {
            case (a1:Array[Float],a2:Array[Float]) => (a1,a2).zipped.map(_+_)
            case (a1:(Any,Any),a2:(Any,Any)) => (cplus[Any](a1._1,a2._1),cplus[Any](a1._2,a2._2))
        }).asInstanceOf[C]
        def cmult[C](collection:C,scalar:Float):C = (collection match {
            case c:Array[Float] => c.map(scalar*_)
            case c:(Any,Any) => (cmult[Any](c._1,scalar),cmult[Any](c._2,scalar))
        }).asInstanceOf[C]
    }

    abstract class LineShape(expir:Array[Float], thirpeaks:Array[(Float,Float)]) {
        import LineShapeHelpers._

        // information of single peak

        // nparams1: number of local parameters that f1 and derivativef1 take
        val nparams1:Int

        // baseline: baseline of the spectrum
        // gparams are global parameters shared with f1
        def baseline(gparams:Array[Float])(x:Float):Float
        // derivative of baseline with respect to gparams
        def dbaseline(gparams:Array[Float])(x:Float):Array[Float]

        // f1: the line shape of a single peak,
        // gparams are global parameters (shared with baseline), params1 are local parameters
        def f1(freq:Float,max:Float)(gparams:Array[Float],params1:Array[Float])(x:Float):Float
        // df1: derivative of f1 with gparams and params1.
        // return_value._1 are derivative with gparams, return_value._2 are with params1
        def df1(freq:Float,max:Float)(gparams:Array[Float],params1:Array[Float])(x:Float):(Array[Float],Array[Float])

        // apply (freq,max) to f1 and df1
        val f1s = thirpeaks.unzip.zipped.map(f1(_,_) _)
        val df1s = thirpeaks.unzip.zipped.map(df1(_,_) _)

        private def splitparams(params:Array[Float]):Array[Array[Float]] = {
            if(nparams1==0)
                thirpeaks.map(j=>Array[Float]())
            else params.sliding(nparams1,nparams1).toArray
        }

        // theoretical spectrum vector
        // formula:
        // thirvec(xi) = sum( j from 0 to m, f1(freqj,maxj)(gparams,params1j)(xi) ) + baseline(gparams)(xi)
        def thirvec(gparams:Array[Float], params:Array[Float]):Array[Float] = {
            // split params into [param1]
            val param1s = splitparams(params)
            // apply (gparams,params1) to f1(freq,max)
            val fparams = (f1s,param1s).zipped.map(_(gparams,_))
            // calculate sum( j from 0 to m, f1(freq,max)(gparams,params1j)(xi) )
            val peaks = xs.map( x => fparams.map(_(x)).reduce(_+_) )
            // calculate baseline(gparams)(xi)
            val blvec = xs.map( baseline(gparams)_ )
            // calculate baseline + peaks
            cplus(peaks,blvec)
        }

        // derivative of theoretical spectrum vector with gparams and params
        // formula:
        // d(thirvec(xi))/d(gparams) = d(baseline(xi))/d(gparams) + sum( j from 0 to m, d(f1(xi))/d(gparams) )
        // d(thirvec(xi))/d(params1j) = d(f1(xi))/d(params1j)
        // def dthirvec(gparams:Array[Float], params:Array[Float]):(Array[Float],Array[Float]) = {
        //     //TODO: implement
        //     val d_baseline = xs.map(dbaseline(gparams)_)
        //     val d_
        // }

        object loss {
            // total loss:
            // formula:
            // euclidean_loss = sum( i from 0 to n, ( thir(xi)-expir(xi) ) ^ 2 )
            def euclidean(gparams:Array[Float], params:Array[Float]):Float = {
                (thirvec(gparams,params),expir).zipped.map((a,b)=>(a-b)*(a-b)).reduce(_+_)
            }

            // d(euclidean)/d(parameters)
            // formula:
            // derivative_gparams = sum( i from 0 to n, 2 * ( thir(xi)-expir(xi) ) * d(thir(xi))/d(gparameters) )
            //                    = sum{ i from 0 to n, 2 * ( thir(xi)-expir(xi) ) * ( derivativeb + sum { j from 0 to m, derivativef1(freqj,maxj)(gparams,param1j)(xi)._1 })}
            // derivative_params  = sum( i from 0 to n, 2 * ( thir(xi)-expir(xi) ) * d(thir(xi))/d(parameters) )
            //                    = [ sum{ i from 0 to n, 2 * ( thir(xi)-expir(xi) ) * derivativef1(freqj,maxj)(gparams,param1j)(xi)._2 } for j from 0 to m ]
            // def deuclidean(gparams:Array[Float], params:Array[Float]):(Array[Float],Array[Float]) = {
            //     // split params into [param1]
            //     val param1s = splitparams(params)
            //     // apply (gparams,params1) to derivativef1(freq,max)
            //     val dparams = (df1s,param1s).zipped.map(_(gparams,_))
            //     // calculate 2 * ( thir(xi)-expir(xi) ) for all xi
            //     val diff2 = (thirvec(gparams,params),expir).zipped.map(2*(_-_))
            //     // calculate the sum{ i from 0 to n,.... } for peaks
            //     def sumi_calculator(d:Float=>(Array[Float],Array[Float])):(Array[Float],Array[Float]) = {
            //         val diff2deriv = (xs,diff2).zipped.map(cmult(d(_),_))
            //         diff2deriv.reduce(cplus[(Array[Float],Array[Float])])
            //     }
            //     val sumi = dparams.map(sumi_calculator)
            //     val dgparams_peaks = sumi.map(_._1).reduce(cplus[Array[Float]])
            //     val dparams_peaks = sumi.flatMap(_._2)
            //     // calculate the sum{ i from 0 to n,.... } for baseline
            //     val dgparams_baseline = (xs,diff2).zipped.map(cmult(baseline(gparams)(_),_))
            //     // final result
            //     ( cplus(dgparams_peaks,dgparams_baseline) , dparams_peaks )
            // }
        }

    }

}
