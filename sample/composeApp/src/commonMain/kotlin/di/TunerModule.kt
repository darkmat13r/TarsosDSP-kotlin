package di

import org.koin.dsl.module
import tuner.Tuner

val tunerModule = module {
    factory <Tuner> {
        Tuner(get())
    }
}