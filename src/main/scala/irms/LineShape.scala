package irms {

    object LineShapeHelpers {
        // vectorization
        // TODO: maybe it's better to set this as a class method of LineShape ?
        def vec[R](f:Float=>R):Seq[R] = X.xs.map(f)

        // calculate elementwise plus/multiplies of two Seq/tuple2:
        // only tuple2, Seq and Float are supported now...
        def cplus[C](a1:C,a2:C):C = ((a1,a2).asInstanceOf[(Any,Any)] match {
            case (a1:Seq[Any],a2:Seq[Any]) => (a1,a2).zipped.map(cplus[Any])
            case (a1:(Any,Any),a2:(Any,Any)) => (cplus[Any](a1._1,a2._1),cplus[Any](a1._2,a2._2))
            case (a1:Float,a2:Float) => a1+a2
        }).asInstanceOf[C]
        def cmult[C](collection:C,scalar:Float):C = (collection match {
            case c:Seq[Any] => c.map(cmult[Any](_,scalar))
            case c:(Any,Any) => (cmult[Any](c._1,scalar),cmult[Any](c._2,scalar))
            case c:Float => scalar*c
        }).asInstanceOf[C]

    }

    abstract class LineShape(peaks:Seq[(Float,Float)]) {

        // nparams1: number of local parameters that f1 and derivativef1 take
        val nparams1:Int

        // baseline: baseline of the spectrum
        // gparams are global parameters shared with f1
        def baseline(gparams:Seq[Float])(x:Float):Float

        // f1: the line shape of a single peak,
        // gparams are global parameters (shared with baseline), params1 are local parameters
        def f1(freq:Float,max:Float,gparams:Seq[Float],params1:Seq[Float])(x:Float):Float

        // evaluate f(x) for each theoretical peak
        // TODO: maybe it's better to set this as private?
        def fxs[R](f:(Float,Float,Seq[Float],Seq[Float])=>Float=>R, gparams:Seq[Float], params:Seq[Float])(x:Float):Seq[R] = {
            val param1s = if(nparams1==0) peaks.map(j=>Seq[Float]())
                          else params.sliding(nparams1,nparams1).toSeq
            val (freqs,maxs) = peaks.unzip
            (freqs,maxs,param1s).zipped.map( (a,b,c)=>f(a,b,gparams,c)(x) )
        }

        // theoretical spectrum vector
        // formula:
        // thir(xi) = sum( j from 0 to m, f1(freqj,maxj)(gparams,params1j)(xi) ) + baseline(gparams)(xi)
        def thir(gparams:Seq[Float], params:Seq[Float])(x:Float):Float = {
            fxs(f1,gparams,params)(x).reduce(_+_) + baseline(gparams)(x)
        }

    }

    trait LineShapeWithLoss extends LineShape {
        // Loss functions
        // The difference between loss1 and loss is, loss is the final loss value that
        // will be used by users, while loss1 is a single frame loss which might be
        // adjusted further by the class
        def loss1(expir:Seq[Float])(gparams:Seq[Float], params:Seq[Float]):Float
        def loss(expir:Seq[Float])(gparams:Seq[Float], params:Seq[Float]):Float = loss1(expir)(gparams,params)
    }

    trait OptimizedLoss extends LineShapeWithLoss {

        // derivative of baseline with respect to gparams
        def dbaseline(gparams:Seq[Float])(x:Float):Seq[Float]

        // df1: derivative of f1 with gparams and params1.
        // return_value._1 are derivative with gparams, return_value._2 are with params1
        def df1(freq:Float,max:Float,gparams:Seq[Float],params1:Seq[Float])(x:Float):(Seq[Float],Seq[Float])

        // dloss1: derivative of loss1 w.r.t. gparams and params
        // return_value._1 are derivative with gparams, return_value._2 are with params
        def dloss1(expir:Seq[Float])(gparams:Seq[Float],params:Seq[Float]):(Seq[Float],Seq[Float])

        // derivative of theoretical spectrum vector with gparams and params
        // formula:
        // d(thir(xi))/d(gparams) = d(baseline(xi))/d(gparams) + sum( j from 0 to m, d(f1(xi))/d(gparams) )
        // d(thir(xi))/d(params1j) = d(f1(xi))/d(params1j)
        def dthir(gparams:Seq[Float], params:Seq[Float])(x:Float):(Seq[Float],Seq[Float]) = {
            import LineShapeHelpers.cplus
            val dbaseline_dgparams = dbaseline(gparams)(x)
            val (df1s_dgparams,df1s_dparams1s) = fxs(df1,gparams,params)(x).unzip
            val dthir_dgparams = cplus(dbaseline_dgparams,df1s_dgparams.reduce(cplus[Seq[Float]]))
            (dthir_dgparams,df1s_dparams1s.flatMap(a=>a))
        }

        // the function that optimize parameters to minimize loss1
        // the return value is optimized (allparams,f)
        def optimize(f:(Seq[Float])=>Float, df:(Seq[Float])=>Seq[Float], initial:Seq[Float]):(Seq[Float],Float)

        // the final loss function after optimization
        override def loss(expir:Seq[Float])(gparams:Seq[Float], params:Seq[Float]):Float = {
            // adapt loss1 and dloss1 to optimize
            val split = gparams.length
            val allparams = gparams++params
            def loss1opt(allparams:Seq[Float]) = (loss1(expir) _).tupled(allparams.splitAt(split))
            def dloss1opt(allparams:Seq[Float]) = {
                val (dgparams,dparams) = (dloss1(expir) _).tupled(allparams.splitAt(split))
                dgparams ++ dparams
            }
            optimize(loss1opt,dloss1opt,allparams)._2
        }
    }

    trait EuclideanLoss extends OptimizedLoss {
        import LineShapeHelpers._

        // formula:
        // sum( i from 0 to n, ( thir(xi)-expir(xi) ) ^ 2 )
        def loss1(expir:Seq[Float])(gparams:Seq[Float], params:Seq[Float]):Float = {
            (vec(thir(gparams,params)),expir).zipped.map((a,b)=>(a-b)*(a-b)).reduce(_+_)
        }

        // d(euclidean)/d(parameters)
        // formula:
        // sum( i from 0 to n, 2 * ( thir(xi)-expir(xi) ) * d(thir(xi))/d(parameters) )
        def dloss1(expir:Seq[Float])(gparams:Seq[Float],params:Seq[Float]):(Seq[Float],Seq[Float]) = {
            def d1(x:Float,expir1:Float) = cmult(dthir(gparams,params)(x), 2*(thir(gparams,params)(x)-expir1))
            (X.xs,expir).zipped.map(d1).reduce(cplus[(Seq[Float], Seq[Float])])
        }

    }

}
