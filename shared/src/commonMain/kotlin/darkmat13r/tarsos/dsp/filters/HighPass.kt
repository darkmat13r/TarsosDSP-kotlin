package darkmat13r.tarsos.dsp.filters

import kotlin.math.PI
import kotlin.math.exp

class HighPass(override val frequency : Float, override val sampleRate : Float, override val overlap : Int) : IIRFilter(frequency, sampleRate, overlap) {
    override fun calcCoeff() {
        val fracFreq = frequency/sampleRate
        val x = exp(-2 * PI * fracFreq)
        aCoeff = floatArrayOf(
            ((1 + x) / 2).toFloat(),
            (-(1 + x) / 2).toFloat()
        )
        bCoeff = floatArrayOf(x.toFloat())
    }
}