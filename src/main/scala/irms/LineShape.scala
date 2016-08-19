package irms {

    object LineShapeHelpers {
        // vectorization
        // TODO: maybe it's better to set this as a class method of LineShape ?
        def vec[R](f:Double=>R):Seq[R] = X.xs.map(f)

        // calculate elementwise plus/multiplies of two Seq/tuple2:
        // only tuple2, Seq and Double are supported now...
        def cplus[C](a1:C,a2:C):C = ((a1,a2).asInstanceOf[(Any,Any)] match {
            case (a1:Seq[Any],a2:Seq[Any]) => (a1,a2).zipped.map(cplus[Any])
            case (a1:(Any,Any),a2:(Any,Any)) => (cplus[Any](a1._1,a2._1),cplus[Any](a1._2,a2._2))
            case (a1:Double,a2:Double) => a1+a2
        }).asInstanceOf[C]
        def cmult[C](collection:C,scalar:Double):C = (collection match {
            case c:Seq[Any] => c.map(cmult[Any](_,scalar))
            case c:(Any,Any) => (cmult[Any](c._1,scalar),cmult[Any](c._2,scalar))
            case c:Double => scalar*c
        }).asInstanceOf[C]

    }

    abstract class LineShape(peaks:Seq[(Double,Double)]) {

        // nparams1: number of local parameters that f1 and derivativef1 take
        val nparams1:Int

        // baseline: baseline of the spectrum
        // gparams are global parameters shared with f1
        def baseline(gparams:Seq[Double])(x:Double):Double

        // f1: the line shape of a single peak,
        // gparams are global parameters (shared with baseline), params1 are local parameters
        def f1(freq:Double,intensity:Double,gparams:Seq[Double],params1:Seq[Double])(x:Double):Double

        // evaluate f(x) for each theoretical peak
        // TODO: maybe it's better to set this as private?
        def fxs[R](f:(Double,Double,Seq[Double],Seq[Double])=>Double=>R, gparams:Seq[Double], params:Seq[Double])(x:Double):Seq[R] = {
            val param1s = if(nparams1==0) peaks.map(j=>Seq[Double]())
                          else params.sliding(nparams1,nparams1).toSeq
            val (freqs,intensitys) = peaks.unzip
            (freqs,intensitys,param1s).zipped.map( (a,b,c)=>f(a,b,gparams,c)(x) )
        }

        // theoretical spectrum vector
        // formula:
        // thir(xi) = sum( j from 0 to m, f1(freqj,intensityj)(gparams,params1j)(xi) ) + baseline(gparams)(xi)
        def thir(gparams:Seq[Double], params:Seq[Double])(x:Double):Double = {
            fxs(f1,gparams,params)(x).reduce(_+_) + baseline(gparams)(x)
        }

        // loss function
        def loss(expir:Seq[Double]):Double
    }

    object Loss {

        trait PlainLoss extends LineShape {
            def gparams:Seq[Double]
            def params:Seq[Double]
            // Loss functions
            // The difference between loss1 and loss is, loss is the final loss value that
            // will be used by users, while loss1 is a single frame loss which might be
            // adjusted further by the class
            def loss1(expir:Seq[Double])(gparams:Seq[Double], params:Seq[Double]):Double
            def loss(expir:Seq[Double]):Double = loss1(expir)(gparams,params)
        }

        trait OptimizedLoss extends PlainLoss {

            var final_gparams:Option[Seq[Double]] = None
            var final_params:Option[Seq[Double]] = None

            // derivative of baseline with respect to gparams
            def dbaseline(gparams:Seq[Double])(x:Double):Seq[Double]

            // df1: derivative of f1 with gparams and params1.
            // return_value._1 are derivative with gparams, return_value._2 are with params1
            def df1(freq:Double,intensity:Double,gparams:Seq[Double],params1:Seq[Double])(x:Double):(Seq[Double],Seq[Double])

            // dloss1: derivative of loss1 w.r.t. gparams and params
            // return_value._1 are derivative with gparams, return_value._2 are with params
            def dloss1(expir:Seq[Double])(gparams:Seq[Double],params:Seq[Double]):(Seq[Double],Seq[Double])

            // derivative of theoretical spectrum vector with gparams and params
            // formula:
            // d(thir(xi))/d(gparams) = d(baseline(xi))/d(gparams) + sum( j from 0 to m, d(f1(xi))/d(gparams) )
            // d(thir(xi))/d(params1j) = d(f1(xi))/d(params1j)
            def dthir(gparams:Seq[Double], params:Seq[Double])(x:Double):(Seq[Double],Seq[Double]) = {
                import LineShapeHelpers.cplus
                val dbaseline_dgparams = dbaseline(gparams)(x)
                val (df1s_dgparams,df1s_dparams1s) = fxs(df1,gparams,params)(x).unzip
                val dthir_dgparams = cplus(dbaseline_dgparams,df1s_dgparams.reduce(cplus[Seq[Double]]))
                (dthir_dgparams,df1s_dparams1s.flatMap(a=>a))
            }

            // the function that optimize parameters to minimize loss1
            // the return value is optimized (allparams,f)
            def optimize(f:(Seq[Double])=>Double, df:(Seq[Double])=>Seq[Double], initial:Seq[Double]):(Seq[Double],Double)

            // the final loss function after optimization
            override def loss(expir:Seq[Double]):Double = {
                // adapt loss1 and dloss1 to optimize
                val mygparams = gparams
                val myparams = params
                val split = mygparams.length
                val allparams = mygparams++myparams
                def loss1opt(allparams:Seq[Double]) = (loss1(expir) _).tupled(allparams.splitAt(split))
                def dloss1opt(allparams:Seq[Double]) = {
                    val (dgparams,dparams) = (dloss1(expir) _).tupled(allparams.splitAt(split))
                    dgparams ++ dparams
                }
                val (final_allparams,final_loss) = optimize(loss1opt,dloss1opt,allparams)
                val final_allparams_splited = final_allparams.splitAt(split)
                final_gparams = Some(final_allparams_splited._1)
                final_params = Some(final_allparams_splited._2)
                final_loss
            }
        }

        trait RepeatedOptimizedLoss extends OptimizedLoss {
            override def loss(expir:Seq[Double]):Double = {
                var remain_trials = 100
                var min_loss = super.loss(expir)
                var min_gparams = final_gparams
                var min_params = final_gparams
                while(remain_trials>0){
                    remain_trials -= 1
                    val new_loss = super.loss(expir)
                    if(new_loss>min_loss) {
                        final_gparams = min_gparams
                        min_params = min_params
                    } else {
                        min_loss = new_loss
                    }
                }
                min_loss
            }
        }

        object Euclidean{
            trait Plain extends PlainLoss {
                import LineShapeHelpers._
                // sum( i from 0 to n, ( thir(xi)-expir(xi) ) ^ 2 )
                def loss1(expir:Seq[Double])(gparams:Seq[Double], params:Seq[Double]):Double = {
                    (vec(thir(gparams,params)),expir).zipped.map((a,b)=>(a-b)*(a-b)).reduce(_+_)
                }
            }
            trait Optmized extends Plain with OptimizedLoss {
                import LineShapeHelpers._
                // sum( i from 0 to n, 2 * ( thir(xi)-expir(xi) ) * d(thir(xi))/d(parameters) )
                def dloss1(expir:Seq[Double])(gparams:Seq[Double],params:Seq[Double]):(Seq[Double],Seq[Double]) = {
                    def d1(x:Double,expir1:Double) = cmult(dthir(gparams,params)(x), 2*(thir(gparams,params)(x)-expir1))
                    (X.xs,expir).zipped.map(d1).reduce(cplus[(Seq[Double], Seq[Double])])
                }
            }
        }
    }


    object Functions {
        object ZeroBaseline {
            trait Plain extends LineShape {
                def baseline(gparams:Seq[Double])(x:Double):Double = 0
            }
            trait Optimized extends Plain with Loss.OptimizedLoss {
                def dbaseline(gparams:Seq[Double])(x:Double):Seq[Double] = Seq.fill(gparams.length)(0.0)
            }
        }
    }

    object Optimizer {
        trait BreezeOptimizer extends Loss.OptimizedLoss {
            import breeze.optimize._
            import breeze.linalg._
            protected type DV = DenseVector[Double]
            protected type DF = DiffFunction[DV]
            protected val optimizer:Minimizer[DV,DF]
            protected def seq2densevec(s:Seq[Double]):DV = new DV(s.toArray)
            protected def densevec2seq(v:DV):Seq[Double] = v.data.toSeq
            def optimize(f:(Seq[Double])=>Double, df:(Seq[Double])=>Seq[Double], initial:Seq[Double]):(Seq[Double],Double) = {
                val difffunc = new DF {
                    def calculate(x:DV) = ( f(densevec2seq(x)) , seq2densevec(df(densevec2seq(x))) )
                }
                val optimum = optimizer.minimize(difffunc,seq2densevec(initial))
                val loss = difffunc(optimum)
                (densevec2seq(optimum),loss)
            }
        }
        trait LBFGSOptimizer extends BreezeOptimizer {
            import breeze.optimize._
            override protected val optimizer = new LBFGS[DV]()
        }
        trait SimpleSGDOptimizer extends BreezeOptimizer {
            import breeze.optimize.StochasticGradientDescent._
            override protected val optimizer = new SimpleSGD[DV](0.05)
        }
    }

}
