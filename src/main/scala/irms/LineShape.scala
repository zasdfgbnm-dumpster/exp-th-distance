package irms {

	abstract class LineShape(peaks:Seq[(Double,Double)]) extends Serializable {

		// baseline: baseline of the spectrum
		def baseline(x:Double):Double

		// f1: the line shape of a single peak
		def f1(freq:Double,intensity:Double)(x:Double):Double

		// allow hook on peak information
		def getpeaks = peaks

		// theoretical spectrum
		// formula:
		// thir(x) = sum(j=0->m,f1(freq_j,intensity_j)(x)) + baseline(x)
		def thir(x:Double):Double = getpeaks.map(j=>f1(j._1,j._2)(x)).sum + baseline(x)

		def thir:Seq[Double] = Params.X.xs.map(thir _)

		// loss function
		def loss(expir:Seq[Double]):Double
	}

	trait UniformScaled extends LineShape {
		val scaling:Double
		override def getpeaks = super.getpeaks.map(j=>(scaling*j._1,j._2))
	}

	trait PiecewiseScaled extends LineShape {
		import scala.collection.immutable.SortedMap
		val scalings:SortedMap[Double,Double] //a list of (upper bound,scaling factor)
		override def getpeaks = super.getpeaks.map(j=>{
			val (freq,intensity) = j
			val (_,scaling) = scalings.from(freq).head
			(freq*scaling,intensity)
		})
	}

	trait CosineLoss extends LineShape {

		// dot product
		def dot(vec1:Seq[Double],vec2:Seq[Double]):Double = (vec1,vec2).zipped.map(_*_).sum

		// cos(v1,v2) = v1.v2/sqrt(v1^2*v2^2)
		// loss = 1-cos(thir,expir)
		def loss(expir:Seq[Double]):Double = 1-dot(thir,expir)/scala.math.sqrt(dot(thir,thir)*dot(expir,expir))
	}

	trait ZeroBaseline extends LineShape {
		def baseline(x:Double):Double = 0
	}

	trait Gaussian extends LineShape {
		val w:Double
		def f1(freq:Double,intensity:Double)(x:Double):Double = {
			val d = x - freq
			val scaled_d = d/w
			(intensity/w) * scala.math.exp(-scaled_d*scaled_d)
		}
	}

	trait Lorentzian extends LineShape {
		val w:Double
		def f1(freq:Double,intensity:Double)(x:Double):Double = {
			val d = x - freq
			intensity * w / ( d * d + w * w )
		}
	}
}
