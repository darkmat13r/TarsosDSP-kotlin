package utils

import androidx.compose.runtime.Composable
import androidx.compose.ui.interop.LocalUIViewController

@Composable
actual fun requestAudioPermission(): () -> Unit {
    val viewController = LocalUIViewController.current
    return {

    }
}