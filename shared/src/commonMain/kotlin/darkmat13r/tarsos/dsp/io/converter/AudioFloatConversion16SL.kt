package darkmat13r.tarsos.dsp.io.converter

import darkmat13r.tarsos.dsp.io.AudioFloatConverter

/***************************************************************************
 *
 * 16 bit signed/unsigned, little/big-endian
 *
 **************************************************************************/

// PCM 16 bit, signed, little-endian
class AudioFloatConversion16SL : AudioFloatConverter() {
    override fun toFloatArray(
        inBuff: ByteArray,
        inOffset: Int,
        outBuff: FloatArray,
        outOffset: Int,
        outLen: Int
    ): FloatArray {
        var ix = inOffset
        val len = outOffset + outLen
        for (ox in outOffset until len) {
            val low = inBuff[ix++].toInt() and 0xFF
            val high = inBuff[ix++].toInt() shl 8
            outBuff[ox] = (low or high).toShort() * (1.0f / 32767.0f)
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
        var ox = outOffset
        val len = inOffset + inLen
        for (ix in inOffset until len) {
            val x = (inBuff[ix] * 32767.0f).toInt()
            outBuff[ox++] = x.toByte()
            outBuff[ox++] = (x shr 8).toByte()
        }
        return outBuff
    }
}