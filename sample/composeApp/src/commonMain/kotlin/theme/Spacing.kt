package ui.theme

import androidx.compose.material3.MaterialTheme
import androidx.compose.runtime.Composable
import androidx.compose.runtime.ReadOnlyComposable
import androidx.compose.runtime.compositionLocalOf
import androidx.compose.ui.unit.dp

class Spacing{
    val micro = 4.dp
    val small = 8.dp
    val normal = 12.dp
    val medium = 16.dp
    val large = 24.dp
    val extraLarge = 36.dp

}
class Dimen {

}

val LocalSpacing  = compositionLocalOf { Spacing() }
val LocalDimen  = compositionLocalOf { Dimen() }

val MaterialTheme.spacing : Spacing
    @Composable
    @ReadOnlyComposable
    get() = LocalSpacing.current


val MaterialTheme.dimen : Dimen
    @Composable
    @ReadOnlyComposable
    get() = LocalDimen.current