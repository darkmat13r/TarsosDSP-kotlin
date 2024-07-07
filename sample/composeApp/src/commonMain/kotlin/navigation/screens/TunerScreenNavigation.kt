package navigation.screens

import androidx.lifecycle.viewmodel.compose.viewModel
import androidx.navigation.NavController
import androidx.navigation.NavGraphBuilder
import androidx.navigation.NavOptions
import androidx.navigation.compose.composable
import ui.home.HomeScreen
import ui.home.TunerViewModel

const val homeScreenRoute = "home"

fun NavController.navigateToHome(navOptions: NavOptions? = null) {
    navigate(homeScreenRoute)
}


fun NavGraphBuilder.homeScreen() {
    composable(homeScreenRoute) {
        val tunerViewModel = viewModel {
            TunerViewModel()
        }
        HomeScreen(viewModel = tunerViewModel)
    }
}