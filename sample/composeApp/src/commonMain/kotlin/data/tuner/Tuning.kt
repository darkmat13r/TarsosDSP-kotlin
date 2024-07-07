package data.tuner

import utils.format

data class Tuning(
    val note: ChromaticScale? = null,
    val frequency: Float = -1f,
    val deviation: TuningDeviationResult = TuningDeviationResult.NotDetected
) {

    val formattedFrequency by lazy { ChromaticScale.FREQUENCY_FORMAT.format(frequency.toDouble()) }

    fun getTone(): String {
        requireNotNull(note)

        return note.tone
    }


}
