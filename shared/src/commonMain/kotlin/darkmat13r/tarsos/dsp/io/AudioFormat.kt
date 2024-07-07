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
package darkmat13r.tarsos.dsp.io

import kotlin.math.abs /*
 * Copyright (c) 1999, 2007, Oracle and/or its affiliates. All rights reserved.
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER.
 *
 * This code is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License version 2 only, as
 * published by the Free Software Foundation.  Oracle designates this
 * particular file as subject to the "Classpath" exception as provided
 * by Oracle in the LICENSE file that accompanied this code.
 *
 * This code is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * version 2 for more details (a copy is included in the LICENSE file that
 * accompanied this code).
 *
 * You should have received a copy of the GNU General Public License version
 * 2 along with this work; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * Please contact Oracle, 500 Oracle Parkway, Redwood Shores, CA 94065 USA
 * or visit www.oracle.com if you need additional information or have any
 * questions.
 */

/**
 * `AudioFormat` is the class that specifies a particular arrangement of data in a sound stream.
 * By examing the information stored in the audio format, you can discover how to interpret the bits in the
 * binary sound data.
 *
 *
 * Every data LineWavelet has an audio format associated with its data stream. The audio format of a source (playback) data LineWavelet indicates
 * what kind of data the data LineWavelet expects to receive for output.  For a target (capture) data LineWavelet, the audio format specifies the kind
 * of the data that can be read from the LineWavelet.
 * Sound files also have audio formats, of course.
 *
 *
 * The `AudioFormat` class accommodates a number of common sound-file encoding techniques, including
 * pulse-code modulation (PCM), mu-law encoding, and a-law encoding.  These encoding techniques are predefined,
 * but service providers can create new encoding types.
 * The encoding that a specific format uses is named by its `encoding` field.
 *
 *
 * In addition to the encoding, the audio format includes other properties that further specify the exact
 * arrangement of the data.
 * These include the number of channels, sample rate, sample size, byte order, frame rate, and frame size.
 * Sounds may have different numbers of audio channels: one for mono, two for stereo.
 * The sample rate measures how many "snapshots" (samples) of the sound pressure are taken per second, per channel.
 * (If the sound is stereo rather than mono, two samples are actually measured at each instant of time: one for the left channel,
 * and another for the right channel; however, the sample rate still measures the number per channel, so the rate is the same
 * regardless of the number of channels.   This is the standard use of the term.)
 * The sample size indicates how many bits are used to store each snapshot; 8 and 16 are typical values.
 * For 16-bit samples (or any other sample size larger than a byte),
 * byte order is important; the bytes in each sample are arranged in
 * either the "little-endian" or "big-endian" style.
 * For encodings like PCM, a frame consists of the set of samples for all channels at a given
 * point in time, and so the size of a frame (in bytes) is always equal to the size of a sample (in bytes) times
 * the number of channels.  However, with some other sorts of encodings a frame can contain
 * a bundle of compressed data for a whole series of samples, as well as additional, non-sample
 * data.  For such encodings, the sample rate and sample size refer to the data after it is decoded into PCM,
 * and so they are completely different from the frame rate and frame size.
 *
 *
 * An `AudioFormat` object can include a set of
 * properties. A property is a pair of key and value: the key
 * is of type `String`, the associated property
 * value is an arbitrary object. Properties specify
 * additional format specifications, like the bit rate for
 * compressed formats. Properties are mainly used as a means
 * to transport additional information of the audio format
 * to and from the service providers. Therefore, properties
 * are ignored in the AudioFormat method.
 *
 *
 * The following table lists some common properties which
 * service providers should use, if applicable:
 *
 *
 * <table border="">
 * <caption>A table with service providers</caption>
 * <tr>
 * <th>Property key</th>
 * <th>Value type</th>
 * <th>Description</th>
</tr> *
 * <tr>
 * <td>&quot;bitrate&quot;</td>
 * <td>[Integer][Integer]</td>
 * <td>average bit rate in bits per second</td>
</tr> *
 * <tr>
 * <td>&quot;vbr&quot;</td>
 * <td>[Boolean][Boolean]</td>
 * <td>`true`, if the file is encoded in variable bit
 * rate (VBR)</td>
</tr> *
 * <tr>
 * <td>&quot;quality&quot;</td>
 * <td>[Integer][Integer]</td>
 * <td>encoding/conversion quality, 1..100</td>
</tr> *
</table> *
 *
 *
 * Vendors of service providers (plugins) are encouraged
 * to seek information about other already established
 * properties in third party plugins, and follow the same
 * conventions.
 *
 * @author Kara Kytle
 * @author Florian Bomers
 * @since 1.3
 */
class AudioFormat
/**
 * Constructs an `AudioFormat` with the given parameters.
 * The encoding specifies the convention used to represent the data.
 * The other parameters are further explained in the
 * @param encoding                  the audio encoding technique
 * @param sampleRate                the number of samples per second
 * @param sampleSizeInBits  the number of bits in each sample
 * @param channels                  the number of channels (1 for mono, 2 for stereo, and so on)
 * @param frameSize                 the number of bytes in each frame
 * @param frameRate                 the number of frames per second
 * @param bigEndian                 indicates whether the data for a single sample
 * is stored in big-endian byte order (`false`
 * means little-endian)
 */(
    /**
     * The audio encoding technique used by this format.
     */
     val encoding: Encoding,
    /**
     * The number of samples played or recorded per second, for sounds that have this format.
     */
    val sampleRate: Float,
    /**
     * The number of bits in each sample of a sound that has this format.
     */
    val sampleSizeInBits: Int,
    /**
     * The number of audio channels in this format (1 for mono, 2 for stereo).
     */
    val channels: Int,
    /**
     * The number of bytes in each frame of a sound that has this format.
     */
    val frameSize: Int,
    /**
     * The number of frames played or recorded per second, for sounds that have this format.
     */
    val frameRate: Float,
    /**
     * Indicates whether the audio data is stored in big-endian or little-endian order.
     */
    val isBigEndian: Boolean
) {
    // INSTANCE VARIABLES

    /**
     * Obtains the sample rate.
     * For compressed formats, the return value is the sample rate of the uncompressed
     * audio data.
     * When this AudioFormat is used for queries capabilities , a sample rate of
     * `AudioSystem.NOT_SPECIFIED` means that any sample rate is
     * acceptable. `AudioSystem.NOT_SPECIFIED` is also returned when
     * the sample rate is not defined for this audio format.
     * @return the number of samples per second,
     * or `AudioSystem.NOT_SPECIFIED`
     *
     * @see .getFrameRate
     */

    /**
     * Obtains the size of a sample.
     * For compressed formats, the return value is the sample size of the
     * uncompressed audio data.
     * When this AudioFormat is used for queries or capabilities , a sample size of
     * `AudioSystem.NOT_SPECIFIED` means that any sample size is
     * acceptable. `AudioSystem.NOT_SPECIFIED` is also returned when
     * the sample size is not defined for this audio format.
     * @return the number of bits in each sample,
     * or `AudioSystem.NOT_SPECIFIED`
     *
     * @see .getFrameSize
     */

    /**
     * Obtains the number of channels.
     * When this AudioFormat is used for queries  or capabilities , a return value of
     * `AudioSystem.NOT_SPECIFIED` means that any (positive) number of channels is
     * acceptable.
     * @return The number of channels (1 for mono, 2 for stereo, etc.),
     * or `AudioSystem.NOT_SPECIFIED`
     */

    /**
     * Obtains the frame size in bytes.
     * When this AudioFormat is used for queries or capabilities, a frame size of
     * `AudioSystem.NOT_SPECIFIED` means that any frame size is
     * acceptable. `AudioSystem.NOT_SPECIFIED` is also returned when
     * the frame size is not defined for this audio format.
     * @return the number of bytes per frame,
     * or `AudioSystem.NOT_SPECIFIED`
     *
     * @see .getSampleSizeInBits
     */

    /**
     * Obtains the frame rate in frames per second.
     * When this AudioFormat is used for queries or capabilities , a frame rate of
     * `AudioSystem.NOT_SPECIFIED` means that any frame rate is
     * acceptable. `AudioSystem.NOT_SPECIFIED` is also returned when
     * the frame rate is not defined for this audio format.
     * @return the number of frames per second,
     * or `AudioSystem.NOT_SPECIFIED`
     *
     * @see .getSampleRate
     */

    /**
     * Indicates whether the audio data is stored in big-endian or little-endian
     * byte order.  If the sample size is not more than one byte, the return value is
     * irrelevant.
     * @return `true` if the data is stored in big-endian byte order,
     * `false` if little-endian
     */


    /** The set of properties  */
    private var properties: HashMap<String, Any> = HashMap()


    /**
     * Constructs an `AudioFormat` with the given parameters.
     * The encoding specifies the convention used to represent the data.
     * The other parameters are further explained in the
     * @param encoding         the audio encoding technique
     * @param sampleRate       the number of samples per second
     * @param sampleSizeInBits the number of bits in each sample
     * @param channels         the number of channels (1 for mono, 2 for
     * stereo, and so on)
     * @param frameSize        the number of bytes in each frame
     * @param frameRate        the number of frames per second
     * @param bigEndian        indicates whether the data for a single sample
     * is stored in big-endian byte order
     * (`false` means little-endian)
     * @param properties       a `Map<String,Object>` object
     * containing format properties
     *
     * @since 1.5
     */
    constructor(
        encoding: Encoding, sampleRate: Float,
        sampleSizeInBits: Int, channels: Int,
        frameSize: Int, frameRate: Float,
        bigEndian: Boolean, properties: Map<String, Any>
    ) : this(
        encoding, sampleRate, sampleSizeInBits, channels,
        frameSize, frameRate, bigEndian
    ) {
        this.properties.putAll(properties)
    }


    /**
     * Constructs an `AudioFormat` with a linear PCM encoding and
     * the given parameters.  The frame size is set to the number of bytes
     * required to contain one sample from each channel, and the frame rate
     * is set to the sample rate.
     *
     * @param sampleRate                the number of samples per second
     * @param sampleSizeInBits  the number of bits in each sample
     * @param channels                  the number of channels (1 for mono, 2 for stereo, and so on)
     * @param signed                    indicates whether the data is signed or unsigned
     * @param bigEndian                 indicates whether the data for a single sample
     * is stored in big-endian byte order (`false`
     * means little-endian)
     */
    constructor(
        sampleRate: Float, sampleSizeInBits: Int,
        channels: Int, signed: Boolean, bigEndian: Boolean
    ) : this(
        (if (signed == true) Encoding.PCM_SIGNED else Encoding.PCM_UNSIGNED),
        sampleRate,
        sampleSizeInBits,
        channels,
        if ((channels == NOT_SPECIFIED || sampleSizeInBits == NOT_SPECIFIED)) NOT_SPECIFIED else ((sampleSizeInBits + 7) / 8) * channels,
        sampleRate,
        bigEndian
    )



    /**
     * Obtain an unmodifiable map of properties.
     * The concept of properties is further explained in
     * the.
     *
     * @return a `Map<String,Object>` object containing
     * all properties. If no properties are recognized, an empty map is
     * returned.
     *
     * @see .getProperty
     * @since 1.5
     */
    fun properties(): Map<String, Any> {
        val ret = properties.toMap()
        return ret
    }


    /**
     * Obtain the property value specified by the key.
     * The concept of properties is further explained in
     * the.
     *
     *
     * If the specified property is not defined for a
     * particular file format, this method returns
     * `null`.
     *
     * @param key the key of the desired property
     * @return the value of the property with the specified key,
     * or `null` if the property does not exist.
     *
     * @see .properties
     * @since 1.5
     */
    fun getProperty(key: String?): Any? {
        if (properties == null) {
            return null
        }
        return properties.get(key)
    }


    /**
     * Indicates whether this format matches the one specified.  To match,
     * two formats must have the same encoding, the same number of channels,
     * and the same number of bits per sample and bytes per frame.
     * The two formats must also have the same sample rate,
     * unless the specified format has the sample rate value `AudioSystem.NOT_SPECIFIED`,
     * which any sample rate will match.  The frame rates must
     * similarly be equal, unless the specified format has the frame rate
     * value `AudioSystem.NOT_SPECIFIED`.  The byte order (big-endian or little-endian)
     * must match if the sample size is greater than one byte.
     *
     * @param format format to test for match
     * @return `true` if this format matches the one specified,
     * `false` otherwise.
     */
    /*
     * $$kk: 04.20.99: i changed the semantics of this.
     */
    fun matches(format: AudioFormat): Boolean {
        if (format.encoding == encoding &&
            ((format.sampleRate == NOT_SPECIFIED.toFloat()) || (format.sampleRate == sampleRate)) &&
            (format.sampleSizeInBits == sampleSizeInBits) &&
            (format.channels == channels &&
                    (format.frameSize == frameSize) &&
                    ((format.frameRate == NOT_SPECIFIED.toFloat()) || (format.frameRate == frameRate)) &&
                    ((format.sampleSizeInBits <= 8) || (format.isBigEndian == isBigEndian)))
        ) return true

        return false
    }


    /**
     * Returns a string that describes the format, such as:
     * "PCM SIGNED 22050 Hz 16 bit mono big-endian".  The contents of the string
     * may vary between implementations of Java Sound.
     *
     * @return a string that describes the format parameters
     */
    override fun toString(): String {
        var sEncoding = ""
        sEncoding = "$encoding "
        val sSampleRate = if (sampleRate == NOT_SPECIFIED.toFloat()) {
            "unknown sample rate, "
        } else {
            "$sampleRate Hz, "
        }
        val sSampleSizeInBits = if (sampleSizeInBits.toFloat() == NOT_SPECIFIED.toFloat()) {
            "unknown bits per sample, "
        } else {
            "" + sampleSizeInBits + " bit, "
        }
        val sChannels = if (channels == 1) {
            "mono, "
        } else if (channels == 2) {
            "stereo, "
        } else {
            if (channels == NOT_SPECIFIED) {
                " unknown number of channels, "
            } else {
                "" + channels + " channels, "
            }
        }
        val sFrameSize = if (frameSize.toFloat() == NOT_SPECIFIED.toFloat()) {
            "unknown frame size, "
        } else {
            "" + frameSize + " bytes/frame, "
        }

        var sFrameRate = ""
        if (abs((sampleRate - frameRate).toDouble()) > 0.00001) {
            sFrameRate =
                if (frameRate == NOT_SPECIFIED.toFloat()) {
                    "unknown frame rate, "
                } else {
                    frameRate.toString() + " frames/second, "
                }
        }

        var sEndian = ""
        if ((encoding == Encoding.PCM_SIGNED || encoding == Encoding.PCM_UNSIGNED)
            && ((sampleSizeInBits > 8)
                    || (sampleSizeInBits == NOT_SPECIFIED))
        ) {
            sEndian = if (isBigEndian) {
                "big-endian"
            } else {
                "little-endian"
            }
        }

        return (sEncoding
                + sSampleRate
                + sSampleSizeInBits
                + sChannels
                + sFrameSize
                + sFrameRate
                + sEndian)
    }

    /**
     * The `Encoding` class  names the  specific type of data representation
     * used for an audio stream.   The encoding includes aspects of the
     * sound format other than the number of channels, sample rate, sample size,
     * frame rate, frame size, and byte order.
     *
     *
     * One ubiquitous type of audio encoding is pulse-code modulation (PCM),
     * which is simply a linear (proportional) representation of the sound
     * waveform.  With PCM, the number stored in each sample is proportional
     * to the instantaneous amplitude of the sound pressure at that point in
     * time.  The numbers are frequently signed or unsigned integers.
     * Besides PCM, other encodings include mu-law and a-law, which are nonlinear
     * mappings of the sound amplitude that are often used for recording speech.
     *
     *
     * You can use a predefined encoding by referring to one of the static
     * objects created by this class, such as PCM_SIGNED or
     * PCM_UNSIGNED.  Service providers can create new encodings, such as
     * compressed audio formats or floating-point PCM samples, and make
     * these available through the `AudioSystem` class.
     *
     *
     * The `Encoding` class is static, so that all
     * `AudioFormat` objects that have the same encoding will refer
     * to the same object (rather than different instances of the same class).
     * This allows matches to be made by checking that two format's encodings
     * are equal.
     *
     * @author Kara Kytle
     * @since 1.3
     */
    class Encoding
    /**
     * Constructs a new encoding.
     * @param name  the name of the new type of encoding
     */(
        /**
         * Encoding name.
         */
        private val name: String
    ) {
        // INSTANCE VARIABLES


        // CONSTRUCTOR


        // METHODS
        /**
         * Finalizes the equals method
         */
        override fun equals(obj: Any?): Boolean {
            if (toString() == null) {
                return (obj != null) && (obj.toString() == null)
            }
            if (obj is Encoding) {
                return toString() == obj.toString()
            }
            return false
        }

        /**
         * Finalizes the hashCode method
         */
        override fun hashCode(): Int {
            if (toString() == null) {
                return 0
            }
            return toString().hashCode()
        }

        /**
         * Provides the `String` representation of the encoding.  This `String` is
         * the same name that was passed to the constructor.  For the predefined encodings, the name
         * is similar to the encoding's variable (field) name.  For example, `PCM_SIGNED.toString()` returns
         * the name "pcm_signed".
         *
         * @return the encoding name
         */
        override fun toString(): String {
            return name
        }

        companion object {
            // ENCODING DEFINES
            /**
             * Specifies signed, linear PCM data.
             */
            val PCM_SIGNED: Encoding = Encoding("PCM_SIGNED")

            /**
             * Specifies unsigned, linear PCM data.
             */
            val PCM_UNSIGNED: Encoding = Encoding("PCM_UNSIGNED")

            /**
             * Specifies u-law encoded data.
             */
            val ULAW: Encoding = Encoding("ULAW")

            /**
             * Specifies a-law encoded data.
             */
            val ALAW: Encoding = Encoding("ALAW")
        }
    } // class Encoding

    companion object {
        const val NOT_SPECIFIED: Int = -1
    }
}
