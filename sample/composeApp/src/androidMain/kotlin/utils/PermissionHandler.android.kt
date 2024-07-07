package utils

import android.util.Log
import androidx.activity.compose.rememberLauncherForActivityResult
import androidx.activity.result.contract.ActivityResultContracts
import androidx.compose.runtime.Composable
import androidx.compose.ui.platform.LocalContext

@Composable
actual fun requestAudioPermission( ): () -> Unit {
    val localContext = LocalContext.current
    val launcher = rememberLauncherForActivityResult(
        ActivityResultContracts.RequestPermission()
    ) { isGranted: Boolean ->
        if (isGranted) {
            // Permission Accepted: Do something
            Log.d("requestAudioPermission","PERMISSION GRANTED")

        } else {
            // Permission Denied: Do something
            Log.d("requestAudioPermission","PERMISSION DENIED")
        }
    }
    return {
        launcher.launch(android.Manifest.permission.RECORD_AUDIO)
    }
}