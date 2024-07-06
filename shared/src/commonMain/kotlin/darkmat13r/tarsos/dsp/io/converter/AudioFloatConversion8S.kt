package darkmat13r.tarsos.dsp.io.converter

import darkmat13r.tarsos.dsp.io.AudioFloatConverter
/***************************************************************************
 *
 * 8 bit signed/unsigned
 *
 **************************************************************************/

// PCM 8 bit, signed
class AudioFloatConversion8S : AudioFloatConverter() {
    override fun toFloatArray(
        inBuff: ByteArray,
        inOffset: Int,
        outBuff: FloatArray,
        outOffset: Int,
        outLen: Int
    ): FloatArray {
        var ix = inOffset
        var ox = outOffset
        for (i in 0 until outLen) {
            outBuff[ox++] = inBuff[ix++].toFloat() * (1.0f / 127.0f)
        }
        return outBuff
    }

    override fun toByteArray(
        inBuff: FloatArray,
        inOffset: Int,
        inLen: Int,
        outBuff: ByteArray,
        outOffset: Int
    ): ByteArray {
        var ix = inOffset
        var ox = outOffset
        for (i in 0 until inLen) {
            outBuff[ox++] = (inBuff[ix++] * 127.0f).toInt().toByte()
        }
        return outBuff
    }
}