package darkmat13r.tarsos.dsp.filters

import kotlin.math.exp
import kotlin.math.max
import kotlin.math.pow

/**
 * Four stage low pass filter.
 */
class LowPassFS(
    override val frequency: Float, override val sampleRate: Float,
    override val overlap: Int
) : IIRFilter(
    //minimum frequency is 60Hz!
    max(60f, frequency), sampleRate, overlap
) {


    override fun calcCoeff() {
        val freqFrac = frequency / sampleRate
        val x = exp(-14 * freqFrac)
        aCoeff = floatArrayOf(
            (1f - x).pow(4.0f)
        )
        bCoeff = floatArrayOf(
            4 * x,
            -6 * x * x,
            4 * x * x * x,
            -x * x * x * x
        )
    }
}