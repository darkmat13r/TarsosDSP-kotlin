package darkmat13r.tarsos.dsp.io.converter

import darkmat13r.tarsos.dsp.io.AudioFloatConverter
import darkmat13r.tarsos.dsp.io.AudioFormat

/***************************************************************************
 *
 * LSB Filter, used filter least significant byte in samples arrays.
 *
 * Is used filter out data in lsb byte when SampleSizeInBits is not
 * dividable by 8.
 *
 **************************************************************************/
class AudioFloatLSBFilter(
    private val converter: AudioFloatConverter,
    private var audioFormat: AudioFormat
) : AudioFloatConverter() {
    private var offset = 0

    private var stepsize = 0

    private var mask: Byte = 0

    private var maskBuffer: ByteArray? = null

    init {
        format = audioFormat
        val bits = audioFormat.sampleSizeInBits
        val bigEndian = audioFormat.isBigEndian
        stepsize = (bits + 7) / 8
        offset = if (bigEndian) (stepsize - 1) else 0
        val lsbBits = bits % 8
        mask = when (lsbBits) {
            0 -> 0x00.toByte()
            1 -> 0x80.toByte()
            2 -> 0xC0.toByte()
            3 -> 0xE0.toByte()
            4 -> 0xF0.toByte()
            5 -> 0xF8.toByte()
            6 -> 0xFC.toByte()
            7 -> 0xFE.toByte()
            else -> 0xFF.toByte()
        }
    }

    override fun toFloatArray(
        inBuff: ByteArray,
        inOffset: Int,
        outBuff: FloatArray,
        outOffset: Int,
        outLen: Int
    ): FloatArray {
        if (maskBuffer == null || maskBuffer!!.size < inBuff.size) {
            maskBuffer = ByteArray(inBuff.size)
        }

        maskBuffer?.let { maskBuffer ->
            inBuff.copyInto(maskBuffer, 0, 0, inBuff.size)
            val inOffsetEnd = outLen * stepsize
            for(i in (inOffset + offset) until inOffsetEnd step stepsize){
                maskBuffer[i] = (maskBuffer[i].toInt() and mask.toInt()).toByte()
            }
            return converter.toFloatArray(maskBuffer, inOffset, outBuff, outOffset, outLen)
        } ?: throw IllegalStateException("Mask buffer is not initialized")
    }

    override fun toByteArray(
        inBuff: FloatArray,
        inOffset: Int,
        inLen: Int,
        outBuff: ByteArray,
        outOffset: Int
    ): ByteArray {
        val result = converter.toByteArray(inBuff, inOffset, inLen, outBuff, outOffset)
        val outOffsetEnd = inLen * stepsize
        for(i in (outOffset + offset) until  outOffsetEnd step stepsize){
            outBuff[i] = (outBuff[i].toInt() and mask.toInt()).toByte()
        }
        return result
    }


}