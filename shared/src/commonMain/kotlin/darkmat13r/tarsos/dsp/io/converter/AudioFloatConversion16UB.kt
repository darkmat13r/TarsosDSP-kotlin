package darkmat13r.tarsos.dsp.io.converter

import darkmat13r.tarsos.dsp.io.AudioFloatConverter

class AudioFloatConversion16UB : AudioFloatConverter() {
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
            val x = ((inBuff[ix++].toInt() and 0xFF) shl 8) or (inBuff[ix++].toInt() and 0xFF)
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
            outBuff[ox++] = (x ushr 8).toByte()
            outBuff[ox++] = x.toByte()
        }
        return outBuff
    }
}


