package data.tuner

import utils.format
import kotlin.math.pow

// Frequencies for the 4th octave
val baseFrequencies = mapOf(
    "C" to 261.63,
    "C#" to 277.18,
    "D" to 293.66,
    "D#" to 311.13,
    "E" to 329.63,
    "F" to 349.23,
    "F#" to 369.99,
    "G" to 392.00,
    "G#" to 415.30,
    "A" to 440.00,
    "A#" to 466.16,
    "B" to 493.88
)

class ChromaticScale(
    val tone: String,
    val octave: Int,
    val frequency: Float,
    val semitone: Boolean = false
) {

    val formattedFrequency by lazy { FREQUENCY_FORMAT.format(frequency.toDouble()) }

    companion object {
        const val FREQUENCY_FORMAT = "%.2f"
        private const val TOTAL_OCTAVES = 10

        val notes by lazy {
            val notes = arrayListOf<ChromaticScale>()
            for (octave in 0 until TOTAL_OCTAVES) {
                val multiplier = 2.0.pow((octave - 4).toDouble())// 4th octave is the base
                for ((note, frequency) in baseFrequencies) {
                    notes.add(
                        ChromaticScale(
                            tone = note,
                            octave = octave,
                            frequency = (frequency * multiplier).toFloat(),
                            semitone = note.contains("#")
                        )
                    )
                }
            }
            notes
        }

        fun getFlatTone(tone: String) = when (tone) {
            Tone.C -> Tone.D
            Tone.D -> Tone.E
            Tone.F -> Tone.G
            Tone.G -> Tone.A
            Tone.A -> Tone.B
            else -> throw IllegalArgumentException("Can't convert $tone to flat")
        }

        fun getSolfegeTone(tone: String) = when (tone) {
            Tone.C -> "Do"
            Tone.D -> "Re"
            Tone.E -> "Mi"
            Tone.F -> "Fa"
            Tone.G -> "Sol"
            Tone.A -> "La"
            Tone.B -> "Si"
            else -> throw IllegalArgumentException("Can't convert $tone to Solfege notation")
        }
    }
}
