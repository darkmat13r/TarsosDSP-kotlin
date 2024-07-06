package darkmat13r.tarsos.dsp.io.converter

import darkmat13r.tarsos.dsp.io.AudioFloatConverter

// PCM 16 bit, signed, big-endian
class AudioFloatConversion16SB : AudioFloatConverter() {
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
            val high = inBuff[ix++].toInt() shl 8
            val low = inBuff[ix++].toInt() and 0xFF
            outBuff[ox++] = (high or low).toShort() * (1.0f / 32767.0f)
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
            val x = (inBuff[ix++] * 32767.0f).toInt()
            outBuff[ox++] = (x shr 8).toByte()
            outBuff[ox++] = x.toByte()
        }
        return outBuff
    }
}