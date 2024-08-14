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
/*
*  Copyright (c) 2007 - 2008 by Damien Di Fede <ddf@compartmental.net>
*
*   This program is free software; you can redistribute it and/or modify
*   it under the terms of the GNU Library General Public License as published
*   by the Free Software Foundation; either version 2 of the License, or
*   (at your option) any later version.
*
*   This program is distributed in the hope that it will be useful,
*   but WITHOUT ANY WARRANTY; without even the implied warranty of
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*   GNU Library General Public License for more details.
*
*   You should have received a copy of the GNU Library General Public
*   License along with this program; if not, write to the Free Software
*   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/
package darkmat13r.tarsos.dsp.util.fft

import darkmat13r.tarsos.dsp.util.TWO_PI
import darkmat13r.tarsos.dsp.util.fft.WindowFunction
import kotlin.jvm.JvmOverloads
import kotlin.math.PI
import kotlin.math.cos

/**
 * A Blackman window function.
 *
 * @author Damien Di Fede
 * @author Corban Brook
 * @see [The Blackman Window](http://en.wikipedia.org/wiki/Window_function.Blackman_windows)
 */
class BlackmanWindow
/** Constructs a Blackman window with a default alpha value of 0.16  */ @JvmOverloads constructor(
    private val alpha: Float = 0.16f
) : WindowFunction() {
    /**
     * Constructs a Blackman window.
     *
     * @param alpha The Blackman alpha parameter
     */

    override fun value(length: Int, index: Int): Float {
        val a0 = (1 - this.alpha) / 2f
        val a1 = 0.5f
        val a2 = this.alpha / 2f

        return a0 - a1 * cos(TWO_PI * index / (length - 1)).toFloat() + a2 * cos(4 * PI * index / (length - 1))
            .toFloat()
    }
}

