package darkmat13r.tarsos.dsp.io.converter

import darkmat13r.tarsos.dsp.io.AudioFloatConverter

class AudioFloatConversion24UB : AudioFloatConverter() {
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
            var x = ((inBuff[ix++].toInt() and 0xFF) shl 16) or
                    ((inBuff[ix++].toInt() and 0xFF) shl 8) or
                    (inBuff[ix++].toInt() and 0xFF)
            x -= 0x7FFFFF
            outBuff[ox++] = x * (1.0f / 0x7FFFFF)
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
            var x = (inBuff[ix++] * 0x7FFFFF).toInt()
            x += 0x7FFFFF
            outBuff[ox++] = (x ushr 16).toByte()
            outBuff[ox++] = (x ushr 8).toByte()
            outBuff[ox++] = x.toByte()
        }
        return outBuff
    }
}


