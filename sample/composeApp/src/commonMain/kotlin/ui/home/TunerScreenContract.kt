package ui.home

import data.tuner.Tuning
import mvi.UIEffect
import mvi.UIEvent
import mvi.UIState

interface TunerScreenContract {

    data class State(
        val tuning: Tuning = Tuning()
    ) : UIState

    sealed interface Event : UIEvent {
        data object Start : Event
        data object Stop : Event
    }

    sealed interface Effect : UIEffect {

    }
}