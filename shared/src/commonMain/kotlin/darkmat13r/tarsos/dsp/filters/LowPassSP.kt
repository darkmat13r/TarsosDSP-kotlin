package darkmat13r.tarsos.dsp.filters

import kotlin.math.PI
import kotlin.math.exp

class LowPassSP (
    override val frequency: Float, override val sampleRate: Float,
    override val overlap: Int
) : IIRFilter(frequency, sampleRate, overlap) {

    override fun calcCoeff() {
        val fracFreq: Float = frequency / sampleRate
        val x = exp(-2 * PI * fracFreq).toFloat()
        aCoeff = floatArrayOf(1 - x)
        bCoeff = floatArrayOf(x)
    }

}