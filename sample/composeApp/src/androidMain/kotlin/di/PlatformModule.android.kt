package di

import androidx.lifecycle.LifecycleOwner
import androidx.lifecycle.lifecycleScope
import darkmat13r.tarsos.example.MainActivity
import dev.icerock.moko.permissions.PermissionsController
import org.koin.core.module.Module
import org.koin.dsl.module

actual fun platformModules(): Module = module{
    single {
        PermissionsController(get())
    }
}