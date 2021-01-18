const progress_history_table_option = {
  item: `<div class="col-lg-3 col-md-4">
      <div class="card mb-3 widget-content" style="height: 100px">
        <div class="widget-content-outer">
          <div class="widget-content-wrapper">
            <button class="btn btn-danger mr-3 id" onclick="remove_history()">
              <i class="fas fa-times"></i>
            </button>
            <a class="link" href="">
              <div class="widget-content-left">
                <div class="widget-heading name"></div>
                <div class="widget-subheading">Last Modified:
                  <span class="time"></span>
                </div>
              </div>
            </a>
            <div class="widget-content-right">
              <div class="status widget-numbers">
                <span class="curr"></span>/<span class="total"></span>
              </div>
            </div>
          </div>
        </div>
      </div>

  </div>`,
  valueNames: ['time', 'curr', 'total', 'name', {
    name: 'status',
    attr: 'data-status'
  }, {
    name: 'link',
    attr: 'href'
  }, {
    name: 'id',
    attr: 'data-id'
  }]
};

const progress_history_table = new List('progress-history-table', progress_history_table_option, {});

$.get("/process/process-history", {
  type: 'processing',
  name: '_all',
}, data => {
  const info = [];
  data.forEach(e => {
    e.link = "/process/process.html?id=" + e.id;
    const time = new Date(e.time);
    e.time = time.toLocaleString();
    info.push(e)
  });
  progress_history_table.clear();
  progress_history_table.add(info)
});

const remove_history = () => {
  let target = $(event.target);
  if (target.prop("tagName") === "I") {
    target = target.parent()
  }
  $.ajax({
    url: '/process/process-history',
    data: {id: target.data("id"), action: 'DELETE'},
    type: 'POST',
    success: data => { target.parent().parent().parent().parent().remove();}
  });
};
