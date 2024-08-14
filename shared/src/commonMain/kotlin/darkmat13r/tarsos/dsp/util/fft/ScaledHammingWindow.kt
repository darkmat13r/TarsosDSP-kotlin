/*
*      _______                       _____   _____ _____  
*     |__   __|                     |  __ \ / ____|  __ \ 
*        | | __ _ _ __ ___  ___  ___| |  | | (___ | |__) |
*        | |/ _` | '__/ __|/ _ \/ __| |  | |\___ \|  ___/ 
*        | | (_| | |  \__ \ (_) \__ \ |__| |____) | |     
*        |_|\__,_|_|  |___/\___/|___/_____/|_____/|_|     
*                                                         
* -------------------------------------------------------------
*
* TarsosDSP is developed by Joren Six at IPEM, University Ghent
*  
* -------------------------------------------------------------
*
*  Info: http://0110.be/tag/TarsosDSP
*  Github: https://github.com/JorenSix/TarsosDSP
*  Releases: http://0110.be/releases/TarsosDSP/
*  
*  TarsosDSP includes modified source code by various authors,
*  for credits and info, see README.
* 
*/
package darkmat13r.tarsos.dsp.util.fft

import darkmat13r.tarsos.dsp.util.TWO_PI
import darkmat13r.tarsos.dsp.util.fft.WindowFunction
import kotlin.math.cos

class ScaledHammingWindow : WindowFunction() {
    override fun value(length: Int, index: Int): Float {
        val scale = 1.0 / length.toDouble() / 0.54
        val factor: Double = TWO_PI / length.toDouble()
        return (scale * (25.0 / 46.0 - 21.0 / 46.0 * cos(factor * index))).toFloat()
    }
}
