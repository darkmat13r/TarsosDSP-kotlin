package darkmat13r.tarsos.dsp.util

import kotlinx.coroutines.*
import kotlin.math.log2
import kotlin.math.pow

object ConcurrencyUtils {
    private val coroutineScope = CoroutineScope(Dispatchers.Default)

    private var THREADS_BEGIN_N_1D_FFT_2THREADS = 8192

    private var THREADS_BEGIN_N_1D_FFT_4THREADS = 65536

    private var THREADS_BEGIN_N_2D = 65536

    private var THREADS_BEGIN_N_3D = 65536

    private var NTHREADS = prevPow2(getNumberOfProcessors())

    /**
     * Returns the number of available processors.
     *
     * @return number of available processors
     */
    fun getNumberOfProcessors(): Int {
        return 1
    }

    /**
     * Returns the current number of threads.
     *
     * @return the current number of threads.
     */
    fun getNumberOfThreads(): Int {
        return NTHREADS
    }

    /**
     * Sets the number of threads. If n is not a power-of-two number, then the
     * number of threads is set to the closest power-of-two number less than n.
     *
     * @param n The number of threads
     */
    fun setNumberOfThreads(n: Int) {
        NTHREADS = prevPow2(n)
    }

    /**
     * Returns the minimal size of 1D data for which two threads are used.
     *
     * @return the minimal size of 1D data for which two threads are used
     */
    fun getThreadsBeginN_1D_FFT_2Threads(): Int {
        return THREADS_BEGIN_N_1D_FFT_2THREADS
    }

    /**
     * Returns the minimal size of 1D data for which four threads are used.
     *
     * @return the minimal size of 1D data for which four threads are used
     */
    fun getThreadsBeginN_1D_FFT_4Threads(): Int {
        return THREADS_BEGIN_N_1D_FFT_4THREADS
    }

    /**
     * Returns the minimal size of 2D data for which threads are used.
     *
     * @return the minimal size of 2D data for which threads are used
     */
    fun getThreadsBeginN_2D(): Int {
        return THREADS_BEGIN_N_2D
    }

    /**
     * Returns the minimal size of 3D data for which threads are used.
     *
     * @return the minimal size of 3D data for which threads are used
     */
    fun getThreadsBeginN_3D(): Int {
        return THREADS_BEGIN_N_3D
    }

    /**
     * Sets the minimal size of 1D data for which two threads are used.
     *
     * @param n the minimal size of 1D data for which two threads are used
     */
    fun setThreadsBeginN_1D_FFT_2Threads(n: Int) {
        THREADS_BEGIN_N_1D_FFT_2THREADS = if (n < 512) 512 else n
    }

    /**
     * Sets the minimal size of 1D data for which four threads are used.
     *
     * @param n the minimal size of 1D data for which four threads are used
     */
    fun setThreadsBeginN_1D_FFT_4Threads(n: Int) {
        THREADS_BEGIN_N_1D_FFT_4THREADS = if (n < 512) 512 else n
    }

    /**
     * Sets the minimal size of 2D data for which threads are used.
     *
     * @param n the minimal size of 2D data for which threads are used
     */
    fun setThreadsBeginN_2D(n: Int) {
        THREADS_BEGIN_N_2D = n
    }

    /**
     * Sets the minimal size of 3D data for which threads are used.
     *
     * @param n the minimal size of 3D data for which threads are used
     */
    fun setThreadsBeginN_3D(n: Int) {
        THREADS_BEGIN_N_3D = n
    }

    /**
     * Resets the minimal size of 1D data for which two and four threads are used.
     */
    fun resetThreadsBeginN_FFT() {
        THREADS_BEGIN_N_1D_FFT_2THREADS = 8192
        THREADS_BEGIN_N_1D_FFT_4THREADS = 65536
    }

    /**
     * Resets the minimal size of 2D and 3D data for which threads are used.
     */
    fun resetThreadsBeginN() {
        THREADS_BEGIN_N_2D = 65536
        THREADS_BEGIN_N_3D = 65536
    }

    /**
     * Returns the closest power-of-two number greater than or equal to x.
     *
     * @param x the number to process
     * @return the closest power-of-two number greater than or equal to x
     */
    fun nextPow2(x: Int): Int {
        require(x >= 1) { "x must be greater or equal to 1" }
        return if ((x and (x - 1)) == 0) {
            x // x is already a power-of-two number
        } else {
            var y = x
            y--
            y = y or (y shr 1)
            y = y or (y shr 2)
            y = y or (y shr 4)
            y = y or (y shr 8)
            y = y or (y shr 16)
            y++
            y
        }
    }

    /**
     * Returns the closest power-of-two number less than or equal to x.
     *
     * @param x the number to process
     * @return the closest power-of-two number less than or equal to x
     */
    fun prevPow2(x: Int): Int {
        require(x >= 1) { "x must be greater or equal to 1" }
        return 2.0.pow(log2(x.toDouble()).toInt().toDouble()).toInt()
    }

    /**
     * Checks if x is a power-of-two number.
     *
     * @param x the number to process
     * @return true if x is a power-of-two number
     */
    fun isPowerOf2(x: Int): Boolean {
        return x > 0 && (x and (x - 1)) == 0
    }

    /**
     * Causes the currently executing coroutine to delay for the specified number of milliseconds.
     *
     * @param millis the number of milliseconds to delay
     */
    suspend fun delayMillis(millis: Long) {
        delay(millis)
    }

    /**
     * Launches a coroutine to execute a suspending task.
     *
     * @param task a suspending function for execution
     */
    fun launchTask(task: suspend () -> Unit): Job {
        return coroutineScope.launch {
            task()
        }
    }

    /**
     * Waits for all jobs to complete computation.
     *
     * @param jobs The jobs which need completion.
     */
    suspend fun waitForCompletion(jobs: List<Job>) {
        jobs.joinAll()
    }
}