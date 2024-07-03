package darkmat13r.tarsos.dsp.filters

import darkmat13r.tarsos.dsp.AudioProcessor

/**
 * An Infinite Impulse Response, or IIR, filter is a filter that uses a set of
 * coefficients and previous filtered values to filter a stream of audio. It is
 * an efficient way to do digital filtering. IIRFilter is a general IIRFilter
 * that simply applies the filter designated by the filter coefficients so that
 * sub-classes only have to dictate what the values of those coefficients are by
 * defining the <code>calcCoeff()</code> function. When filling the coefficient
 * arrays, be aware that <code>b[0]</code> corresponds to
 * <code>b<sub>1</sub></code>.
 *
 * @author Damien Di Fede
 *
 */
abstract class IIRFilter(val frequency: Float,
                         val sample: Float,
                         private val overlap: Int) :
    AudioProcessor {

    protected lateinit var aCoeff: FloatArray

    protected lateinit var bCoeff: FloatArray

    protected lateinit var input: FloatArray

    protected lateinit var prevOutput: FloatArray


    init {
        this.calcCoeff()

        check(::aCoeff.isInitialized)
        check(::bCoeff.isInitialized)

        input = FloatArray(aCoeff.size)
        prevOutput = FloatArray(bCoeff.size)
    }

    protected abstract fun calcCoeff()

    override fun processFull(audioFloatBuffer: FloatArray, audioByteBuffer: ByteArray): Boolean {
        process(0, audioFloatBuffer)
        return true
    }

    override fun processOverlapping(
        audioFloatBuffer: FloatArray,
        audioByteBuffer: ByteArray
    ): Boolean {
        process(overlap, audioFloatBuffer)
        return true
    }


    private fun process(offset: Int, audioFloatBuffer: FloatArray) {
        for (i in offset until audioFloatBuffer.size){
            //shift the in array
            input.copyInto(input, 1, 0, input.size - 1)
            input[0] = audioFloatBuffer[i]

            //calculate y based on a and b coefficients
            //and in and out.
            var y = 0f
            for (j in aCoeff.indices){
                y += aCoeff[j] * input[j];
            }
            for (j in bCoeff.indices){
                y += bCoeff[j] * prevOutput[j];
            }

            //shift the out array
            prevOutput.copyInto(prevOutput, 1, 0 , prevOutput.size - 1)
            prevOutput[0] = y

            audioFloatBuffer[i] = y
        }
    }

}