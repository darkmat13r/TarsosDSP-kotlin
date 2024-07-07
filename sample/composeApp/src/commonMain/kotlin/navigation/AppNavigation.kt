package navigation

import androidx.compose.runtime.Composable
import androidx.navigation.compose.NavHost
import androidx.navigation.compose.rememberNavController
import navigation.screens.homeScreen
import navigation.screens.homeScreenRoute

@Composable
fun AppNavigation() {
    val navController = rememberNavController()
    NavHost(
        navController,
        startDestination = homeScreenRoute
    ) {
        homeScreen()
    }
}