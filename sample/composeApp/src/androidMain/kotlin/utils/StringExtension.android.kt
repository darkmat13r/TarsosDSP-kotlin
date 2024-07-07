package utils

actual fun String.format(value: Double): String {
    return String.format(this, value)
}