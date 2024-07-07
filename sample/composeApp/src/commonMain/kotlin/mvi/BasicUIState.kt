package mvi

sealed interface BasicUIState<out T> {
    data class Success<T>(val data: T) : BasicUIState<T>
    data class Progress<T>(val data: Int) : BasicUIState<T>
    data class Error(val message: String? = null) : BasicUIState<Nothing>
    data object Unauthorized : BasicUIState<Nothing>
    data object Loading : BasicUIState<Nothing>
    data object Empty : BasicUIState<Nothing>
    data object Idle : BasicUIState<Nothing>
}


