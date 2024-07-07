package ui.home

import androidx.compose.foundation.layout.Arrangement
import androidx.compose.foundation.layout.Column
import androidx.compose.foundation.layout.fillMaxWidth
import androidx.compose.foundation.layout.padding
import androidx.compose.material3.MaterialTheme
import androidx.compose.material3.Scaffold
import androidx.compose.material3.Text
import androidx.compose.runtime.Composable
import androidx.compose.runtime.LaunchedEffect
import androidx.compose.runtime.collectAsState
import androidx.compose.runtime.getValue
import androidx.compose.runtime.rememberCoroutineScope
import androidx.compose.ui.Alignment
import androidx.compose.ui.Modifier
import androidx.compose.ui.platform.LocalLifecycleOwner
import androidx.lifecycle.Lifecycle
import androidx.lifecycle.LifecycleOwner

import kotlinx.coroutines.launch

import ui.theme.spacing
import utils.requestAudioPermission

@Composable
fun HomeScreen(
    viewModel: TunerViewModel,
) {
    val coroutineScope = rememberCoroutineScope()
    val permissionHandler = requestAudioPermission()
    val state by viewModel.uiState.collectAsState()
    val lifecycle = LocalLifecycleOwner.current.lifecycle

    LaunchedEffect(lifecycle) {
        when (lifecycle.currentState) {
            Lifecycle.State.DESTROYED -> {
                viewModel.setEvent(TunerScreenContract.Event.Stop)
            }

            Lifecycle.State.STARTED -> {
                viewModel.setEvent(TunerScreenContract.Event.Start)
            }

            else -> {}
        }
    }

    LaunchedEffect(true) {
        coroutineScope.launch {
            permissionHandler.invoke()
        }
    }

    Scaffold {
        Column(
            modifier = Modifier.padding(it).fillMaxWidth().padding(MaterialTheme.spacing.medium),
            horizontalAlignment = Alignment.CenterHorizontally,
            verticalArrangement = Arrangement.Center
        ) {
            Text("Pitch Detected", style = MaterialTheme.typography.labelSmall)
            Text("${state.tuning.note}", style = MaterialTheme.typography.headlineMedium)
        }
    }
}