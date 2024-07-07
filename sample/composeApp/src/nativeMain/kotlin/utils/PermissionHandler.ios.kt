package utils

import androidx.compose.runtime.Composable
import androidx.compose.ui.interop.LocalUIViewController
import co.touchlab.kermit.Logger
import platform.AVFAudio.AVAudioSession

@Composable
actual fun requestAudioPermission(): () -> Unit {
    return {
        AVAudioSession.sharedInstance().requestRecordPermission { granted ->
            Logger.i("Microphone permission granted: $granted")
        }
    }
}