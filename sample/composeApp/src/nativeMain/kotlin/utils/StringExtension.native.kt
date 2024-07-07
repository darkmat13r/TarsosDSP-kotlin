package utils

import platform.Foundation.NSNumberFormatter
import platform.Foundation.NSString
import platform.Foundation.stringWithFormat

actual fun String.format(value: Double): String {
    val formatString = this
    return NSString.stringWithFormat(formatString, value) as String
}