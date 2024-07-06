package darkmat13r.tarsos.dsp.io.converter

import darkmat13r.tarsos.dsp.io.AudioFloatConverter

class AudioFloatConversion16UL : AudioFloatConverter() {
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
            val low = inBuff[ix++].toInt() and 0xFF
            val high = inBuff[ix++].toInt() shl 8
            val x = low or high
            outBuff[ox++] = (x - 32767) * (1.0f / 32767.0f)
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
            val x = 32767 + (inBuff[ix++] * 32767.0f).toInt()
            outBuff[ox++] = x.toByte()
            outBuff[ox++] = (x ushr 8).toByte()
        }
        return outBuff
    }
}