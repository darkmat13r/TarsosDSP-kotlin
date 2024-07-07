package darkmat13r.tarsos.dsp.io.converter

import darkmat13r.tarsos.dsp.io.AudioFloatConverter

// PCM 32 bit float, big-endian
class AudioFloatConversion32B : AudioFloatConverter(){
    override fun toFloatArray(
        inBuff: ByteArray,
        inOffset: Int,
        outBuff: FloatArray,
        outOffset: Int,
        outLen: Int
    ): FloatArray {
        for (i in 0 until outLen) {
            val index = inOffset + i * 4
            outBuff[outOffset + i] = bytesToFloat(inBuff, index)
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
        for (i in 0 until inLen) {
            val index = outOffset + i * 4
            floatToBytes(inBuff[inOffset + i], outBuff, index)
        }
        return outBuff
    }

    // Convert 4 bytes to a float in big-endian order
    private fun bytesToFloat(bytes: ByteArray, index: Int): Float {
        val intBits = ((bytes[index].toInt() and 0xff) shl 24) or
                ((bytes[index + 1].toInt() and 0xff) shl 16) or
                ((bytes[index + 2].toInt() and 0xff) shl 8) or
                (bytes[index + 3].toInt() and 0xff)
        return Float.fromBits(intBits)
    }

    // Convert a float to 4 bytes in big-endian order
    private fun floatToBytes(value: Float, bytes: ByteArray, index: Int) {
        val intBits = value.toBits()
        bytes[index] = ((intBits shr 24) and 0xff).toByte()
        bytes[index + 1] = ((intBits shr 16) and 0xff).toByte()
        bytes[index + 2] = ((intBits shr 8) and 0xff).toByte()
        bytes[index + 3] = (intBits and 0xff).toByte()
    }
}