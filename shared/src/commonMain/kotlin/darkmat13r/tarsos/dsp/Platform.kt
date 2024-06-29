package darkmat13r.tarsos.dsp

interface Platform {
    val name: String
}

expect fun getPlatform(): Platform