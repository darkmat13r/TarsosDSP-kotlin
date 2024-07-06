package darkmat13r.tarsos.dsp.io

import darkmat13r.tarsos.dsp.io.converter.AudioFloatConversion16SB
import darkmat13r.tarsos.dsp.io.converter.AudioFloatConversion16SL
import darkmat13r.tarsos.dsp.io.converter.AudioFloatConversion16UB
import darkmat13r.tarsos.dsp.io.converter.AudioFloatConversion16UL
import darkmat13r.tarsos.dsp.io.converter.AudioFloatConversion24SB
import darkmat13r.tarsos.dsp.io.converter.AudioFloatConversion24SL
import darkmat13r.tarsos.dsp.io.converter.AudioFloatConversion24UB
import darkmat13r.tarsos.dsp.io.converter.AudioFloatConversion24UL
import darkmat13r.tarsos.dsp.io.converter.AudioFloatConversion32B
import darkmat13r.tarsos.dsp.io.converter.AudioFloatConversion32L
import darkmat13r.tarsos.dsp.io.converter.AudioFloatConversion32SB
import darkmat13r.tarsos.dsp.io.converter.AudioFloatConversion32SL
import darkmat13r.tarsos.dsp.io.converter.AudioFloatConversion32UB
import darkmat13r.tarsos.dsp.io.converter.AudioFloatConversion32UL
import darkmat13r.tarsos.dsp.io.converter.AudioFloatConversion32xSB
import darkmat13r.tarsos.dsp.io.converter.AudioFloatConversion32xSL
import darkmat13r.tarsos.dsp.io.converter.AudioFloatConversion32xUB
import darkmat13r.tarsos.dsp.io.converter.AudioFloatConversion32xUL
import darkmat13r.tarsos.dsp.io.converter.AudioFloatConversion64B
import darkmat13r.tarsos.dsp.io.converter.AudioFloatConversion64L
import darkmat13r.tarsos.dsp.io.converter.AudioFloatConversion8S
import darkmat13r.tarsos.dsp.io.converter.AudioFloatConversion8U
import darkmat13r.tarsos.dsp.io.converter.AudioFloatLSBFilter

abstract class AudioFloatConverter() {
    open var format: AudioFormat? = null

    companion object {
        private val PCM_FLOAT = AudioFormat.Encoding("PCM_FLOAT")
    }


    abstract fun toFloatArray(
        inBuff: ByteArray,
        inOffset: Int = 0,
        outBuff: FloatArray,
        outOffset: Int = 0,
        outLen: Int = outBuff.size
    ): FloatArray


    abstract fun toByteArray(
        inBuff: FloatArray,
        inOffset: Int = 0,
        inLen: Int = inBuff.size,
        outBuff: ByteArray,
        outOffset: Int = 0,
    ): ByteArray


    object Factory {

        fun create(format: AudioFormat): AudioFloatConverter? {
            var conv: AudioFloatConverter? = null
            if (format.frameSize == 0) return null
            if (format.frameSize != ((format.sampleSizeInBits + 7) / 8) * format.channels) return null

            if (format.encoding == AudioFormat.Encoding.PCM_SIGNED) {
                conv = if (format.isBigEndian) {
                    when {
                        format.sampleSizeInBits <= 8 -> AudioFloatConversion8S()
                        format.sampleSizeInBits <= 16 -> AudioFloatConversion16SB()
                        format.sampleSizeInBits <= 24 -> AudioFloatConversion24SB()
                        format.sampleSizeInBits <= 32 -> AudioFloatConversion32SB()
                        else -> AudioFloatConversion32xSB(((format.sampleSizeInBits + 7) / 8) - 4)
                    }
                } else {
                    when {
                        format.sampleSizeInBits <= 8 -> AudioFloatConversion8S()
                        format.sampleSizeInBits <= 16 -> AudioFloatConversion16SL()
                        format.sampleSizeInBits <= 24 -> AudioFloatConversion24SL()
                        format.sampleSizeInBits <= 32 -> AudioFloatConversion32SL()
                        else -> AudioFloatConversion32xSL(((format.sampleSizeInBits + 7) / 8) - 4)
                    }
                }
            } else if (format.encoding == AudioFormat.Encoding.PCM_UNSIGNED) {
                conv = if (format.isBigEndian) {
                    when {
                        format.sampleSizeInBits <= 8 -> AudioFloatConversion8U()
                        format.sampleSizeInBits <= 16 -> AudioFloatConversion16UB()
                        format.sampleSizeInBits <= 24 -> AudioFloatConversion24UB()
                        format.sampleSizeInBits <= 32 -> AudioFloatConversion32UB()
                        else -> AudioFloatConversion32xUB(((format.sampleSizeInBits + 7) / 8) - 4)
                    }
                } else {
                    when {
                        format.sampleSizeInBits <= 8 -> AudioFloatConversion8U()
                        format.sampleSizeInBits <= 16 -> AudioFloatConversion16UL()
                        format.sampleSizeInBits <= 24 -> AudioFloatConversion24UL()
                        format.sampleSizeInBits <= 32 -> AudioFloatConversion32UL()
                        else -> AudioFloatConversion32xUL(((format.sampleSizeInBits + 7) / 8) - 4)
                    }
                }
            } else if (format.encoding == PCM_FLOAT) {
                conv = when (format.sampleSizeInBits) {
                    32 -> if (format.isBigEndian) AudioFloatConversion32B() else AudioFloatConversion32L()
                    64 -> if (format.isBigEndian) AudioFloatConversion64B() else AudioFloatConversion64L()
                    else -> null
                }
            }

            if ((format.encoding == AudioFormat.Encoding.PCM_SIGNED || format.encoding == AudioFormat.Encoding.PCM_UNSIGNED) &&
                (format.sampleSizeInBits % 8 != 0)
            ) {
                if (conv != null) {
                    conv = AudioFloatLSBFilter(conv, format)
                }
            }

            conv?.format = format
            return conv
        }
    }

}