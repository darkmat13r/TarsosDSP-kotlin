package darkmat13r.tarsos.dsp.io

actual class NativeAudioInputStream(
    override var format: AudioFormat,
    override var frameLength: Long
) : AudioInputStream {
    override fun skip(bytesToSkip: Long): Long {
        TODO("Not yet implemented")
    }

    override fun read(b: ByteArray, off: Int, len: Int): Int {
        TODO("Not yet implemented")
    }

    override fun close() {
        TODO("Not yet implemented")
    }
}