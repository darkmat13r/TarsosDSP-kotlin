package coroutines

import kotlinx.coroutines.CoroutineDispatcher
import kotlinx.coroutines.Dispatchers
import kotlinx.coroutines.Runnable
import platform.Foundation.NSRunLoop
import platform.Foundation.performBlock
import platform.darwin.dispatch_get_global_queue
import platform.posix.QOS_CLASS_BACKGROUND
import kotlin.coroutines.CoroutineContext


private class NsQueueDispatcher : CoroutineDispatcher() {
    override fun dispatch(context: CoroutineContext, block: Runnable) {
        NSRunLoop.mainRunLoop().performBlock { block.run() }
    }
}

actual val IODispatcher: CoroutineDispatcher
    get() = NsQueueDispatcher()