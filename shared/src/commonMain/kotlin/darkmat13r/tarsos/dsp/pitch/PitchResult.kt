package darkmat13r.tarsos.dsp.pitch

sealed interface PitchResult{
     class Pitch : PitchResult{
         /**
          * pitch in Hz
          */
         var hz : Float = 0f
         var probability : Float = -1f
         var pitched : Boolean = false

         override fun toString(): String {
             return "Hz : $hz, Probability $probability,  Pitched $pitched"
         }
     }
    data object None : PitchResult
}