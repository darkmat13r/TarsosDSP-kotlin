package darkmat13r.tarsos.dsp.io.converter

import darkmat13r.tarsos.dsp.io.AudioFloatConverter

// PCM 8 bit, unsigned
class AudioFloatConversion8U : AudioFloatConverter() {
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
            outBuff[ox++] = ((inBuff[ix++].toInt() and 0xFF) - 127) * (1.0f / 127.0f)
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
            outBuff[ox++] = (127 + (inBuff[ix++] * 127.0f).toInt()).toByte()
        }
        return outBuff
    }
}