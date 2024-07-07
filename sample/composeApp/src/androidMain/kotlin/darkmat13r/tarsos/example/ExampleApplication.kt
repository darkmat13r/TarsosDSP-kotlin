package darkmat13r.tarsos.example

import android.app.Application
import di.initKoin
import org.koin.android.ext.koin.androidContext

class ExampleApplication : Application() {

    override fun onCreate() {
        super.onCreate()
        initKoin{
            androidContext(this@ExampleApplication)
        }
    }
}