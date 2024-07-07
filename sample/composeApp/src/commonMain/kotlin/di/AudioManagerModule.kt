package di

import audio.AudioManager
import org.koin.dsl.module

val audioManager = module {
    single {
        AudioManager()
    }
}