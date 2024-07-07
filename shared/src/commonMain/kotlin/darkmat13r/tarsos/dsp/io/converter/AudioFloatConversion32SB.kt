package darkmat13r.tarsos.dsp.io.converter

import darkmat13r.tarsos.dsp.io.AudioFloatConverter

class AudioFloatConversion32SB : AudioFloatConverter() {
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
            val x = ((inBuff[ix++].toInt() and 0xFF) shl 24) or
                    ((inBuff[ix++].toInt() and 0xFF) shl 16) or
                    ((inBuff[ix++].toInt() and 0xFF) shl 8) or
                    (inBuff[ix++].toInt() and 0xFF)
            outBuff[ox++] = x * (1.0f / 0x7FFFFFFF)
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
            val x = (inBuff[ix++] * 0x7FFFFFFF).toInt()
            outBuff[ox++] = (x ushr 24).toByte()
            outBuff[ox++] = (x ushr 16).toByte()
            outBuff[ox++] = (x ushr 8).toByte()
            outBuff[ox++] = x.toByte()
        }
        return outBuff
    }
}

