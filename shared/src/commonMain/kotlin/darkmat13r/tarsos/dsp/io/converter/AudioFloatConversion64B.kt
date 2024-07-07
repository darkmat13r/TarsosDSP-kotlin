package darkmat13r.tarsos.dsp.io.converter

import darkmat13r.tarsos.dsp.io.AudioFloatConverter

// PCM 64 bit float, big-endian
class AudioFloatConversion64B : AudioFloatConverter() {
    private var doubleBuff: DoubleArray? = null

    override fun toFloatArray(
        inBuff: ByteArray,
        inOffset: Int,
        outBuff: FloatArray,
        outOffset: Int,
        outLen: Int
    ): FloatArray {
        if (doubleBuff == null || doubleBuff!!.size < outLen + outOffset) {
            doubleBuff = DoubleArray(outLen + outOffset)
        }
        for (i in 0 until outLen) {
            val index = inOffset + i * 8
            doubleBuff!![outOffset + i] = bytesToDouble(inBuff, index)
        }
        for (i in outOffset until outOffset + outLen) {
            outBuff[i] = doubleBuff!![i].toFloat()
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
        if (doubleBuff == null || doubleBuff!!.size < inOffset + inLen) {
            doubleBuff = DoubleArray(inOffset + inLen)
        }
        for (i in 0 until inLen) {
            doubleBuff!![inOffset + i] = inBuff[i].toDouble()
            val index = outOffset + i * 8
            doubleToBytes(doubleBuff!![inOffset + i], outBuff, index)
        }
        return outBuff
    }

    // Convert 8 bytes to a double in big-endian order
    private fun bytesToDouble(bytes: ByteArray, index: Int): Double {
        val longBits = ((bytes[index].toLong() and 0xff) shl 56) or
                ((bytes[index + 1].toLong() and 0xff) shl 48) or
                ((bytes[index + 2].toLong() and 0xff) shl 40) or
                ((bytes[index + 3].toLong() and 0xff) shl 32) or
                ((bytes[index + 4].toLong() and 0xff) shl 24) or
                ((bytes[index + 5].toLong() and 0xff) shl 16) or
                ((bytes[index + 6].toLong() and 0xff) shl 8) or
                (bytes[index + 7].toLong() and 0xff)
        return Double.fromBits(longBits)
    }

    // Convert a double to 8 bytes in big-endian order
    private fun doubleToBytes(value: Double, bytes: ByteArray, index: Int) {
        val longBits = value.toBits()
        bytes[index] = ((longBits shr 56) and 0xff).toByte()
        bytes[index + 1] = ((longBits shr 48) and 0xff).toByte()
        bytes[index + 2] = ((longBits shr 40) and 0xff).toByte()
        bytes[index + 3] = ((longBits shr 32) and 0xff).toByte()
        bytes[index + 4] = ((longBits shr 24) and 0xff).toByte()
        bytes[index + 5] = ((longBits shr 16) and 0xff).toByte()
        bytes[index + 6] = ((longBits shr 8) and 0xff).toByte()
        bytes[index + 7] = (longBits and 0xff).toByte()
    }
}