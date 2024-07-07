import androidx.compose.runtime.Composable
import navigation.AppNavigation
import org.jetbrains.compose.ui.tooling.preview.Preview
import theme.AppTheme

@Composable
@Preview
fun App() {
    AppTheme {
        AppNavigation()
    }
}