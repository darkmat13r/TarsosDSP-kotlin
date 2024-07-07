package ui.home

import androidx.lifecycle.viewModelScope
import coroutines.IODispatcher
import kotlinx.coroutines.flow.launchIn
import kotlinx.coroutines.flow.onEach
import kotlinx.coroutines.launch
import mvi.BaseViewModel
import org.koin.core.component.KoinComponent
import org.koin.core.component.inject
import tuner.Tuner

class TunerViewModel :
    BaseViewModel<TunerScreenContract.Event, TunerScreenContract.State, TunerScreenContract.Effect>(),
    KoinComponent {

    private val tuner: Tuner by inject<Tuner>()

    init {
        viewModelScope.launch {
            setupState()
        }
    }

    private fun setupState() {
        tuner.state
            .onEach { tuning ->
                setState {
                    copy(tuning = tuning)
                }
            }.launchIn(viewModelScope)
    }

    override fun createInitialState(): TunerScreenContract.State {
        return TunerScreenContract.State()
    }


    override fun handleEvent(event: TunerScreenContract.Event) {
        when (event) {
            TunerScreenContract.Event.Start -> {
                viewModelScope.launch(IODispatcher) {
                    tuner.startRecording()
                }
            }

            TunerScreenContract.Event.Stop -> viewModelScope.launch(IODispatcher) {
                tuner.stopRecording()
            }
        }
    }
}